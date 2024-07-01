use clap::{command, Arg};
use core::fmt;
use std::io::Read;
use std::io::{self};
use std::io::BufRead;

#[derive(Debug)]
struct Probe {
    name: String,
    sd: f64,
}

#[derive(Debug)]
struct TopProbes {
    top_probes: Vec<Probe>,
    n_probes_failed_qc: usize,
    n_total_probes: usize,
    n_samples: usize,
}

impl fmt::Display for Probe {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\t{}", self.name, self.sd)
    }
}

#[derive(Debug)]
enum InputSource {
    Stdin,
    File(String),
}

fn main() {

    // Argument Parsing
    let matches = command!() // requires `cargo` feature
        .subcommand(
            clap::Command::new("identify")
            .about("Identify the most variable probes in a methylation dataset")
            .arg(Arg::new("methdata")
                .help("A gzipped tsv with methylation data where first column has probe names, all other columns are samples. Values represent beta values")
                .required(false)
            )
            .arg(
                Arg::new("delim")
                .short('d')
                .long("delim")
                .default_value("\t")
                .value_parser(clap::value_parser!(String))
                .help("File delimiter")
            )
            .arg(
                Arg::new("nprobes")
                .short('n')
                .long("nprobes")
                .default_value("10")
                .value_parser(clap::value_parser!(usize))
                .help("How many probes should be returned")
            )
            .arg(
                Arg::new("proportion")
                .short('c')
                .long("coverage")
                .default_value("0.02")
                .value_parser(clap::value_parser!(f64))
                .help("Maximum permissible proportion of samples with non-null values before probe fails QC")
            )
        )
        .subcommand(
            clap::Command::new("select")
            .about("Subset a methylation array dataset by sample IDs")
            .arg(Arg::new("methdata_select")
                .help("A gzipped tsv with methylation data where first column has probe names, all other columns are samples. Must have header line")
                .required(false)
            )
            .arg(
                Arg::new("sample_file")
                .short('s')
                .long("samples")
                .help("A plain text file with 1 sample identifier per line. No header")
                .required(true)
            )
            .arg(
                Arg::new("delim_select")
                .short('d')
                .long("delim")
                .default_value("\t")
                .value_parser(clap::value_parser!(String))
                .help("File delimiter")
            )
            
        )
        .about("Pick the most variable probes from a methylation dataset")
        .get_matches();



    // subcommand: SELECT
    if let Some(matches) = matches.subcommand_matches("select") { 
        let delim = matches
        .get_one::<String>("delim_select")
        .expect("Valid file delimiter")
        .as_bytes()[0];

        let methdata_path_option = matches.get_one::<String>("methdata_select");
        let sample_path = matches.get_one::<String>("sample_file").expect("Sample Path");

        let reader: csv::Reader<Box<dyn Read>> = get_reader_from_stdin_or_file( methdata_path_option, delim);
        let samples: Vec<String> = get_selected_samples(sample_path.to_string());
        
        print_selected_samples(reader, samples)
     }    
     // subcommand: IDENTIFY 
     else if let Some(matches) = matches.subcommand_matches("identify") { 
        let delim = matches
            .get_one::<String>("delim")
            .expect("Valid file delimiter")
            .as_bytes()[0];

        let methdata_path_option = matches.get_one::<String>("methdata");

        let reader = get_reader_from_stdin_or_file( methdata_path_option, delim);
        
        // Config
        let n_probes_to_return: usize = *matches
            .get_one::<usize>("nprobes")
            .expect("Failed to parse nprobes argument");

        let max_proportion_missing_samples = *matches
            .get_one::<f64>("proportion")
            .expect("Failed to parse coverage argument");

        let output_probenames_only = true;

        // Identify Most Variable Probes
        let results =
            get_most_variable_probes(reader, n_probes_to_return, max_proportion_missing_samples);

        // Print Results
        summarise_identify_results(results, max_proportion_missing_samples, n_probes_to_return, output_probenames_only)
    }  
    else {
        panic!("Please supply a valid subcommand")
    }

}

fn summarise_identify_results(results:TopProbes, max_proportion_missing_samples: f64, n_probes_to_return: usize, output_probenames_only: bool){
     // Print Out Debug Information
     eprintln!(
        "Looked through {} probes and {} samples",
        results.n_total_probes, results.n_samples
    );
    eprintln!(
        "{} / {} ({:.0}%) of probes passed QC (>= {:.1}% sample coverage):",
        results.n_total_probes - results.n_probes_failed_qc,
        results.n_total_probes,
        (results.n_total_probes - results.n_probes_failed_qc) as f64
            / results.n_total_probes as f64
            * 100.0,
        100.0 - max_proportion_missing_samples * 100.0,
    );
    if results.top_probes.len() < n_probes_to_return {
        eprintln!(
            "Failed to find the desired {} probes to return. This is usually because too many probes fail QC",
            n_probes_to_return
        );
    }
    eprintln!(
        "Identified the {} most variable probes",
        results.top_probes.len()
    );
    eprint_hr();
    print_probe_vector(&results.top_probes, output_probenames_only);
}

fn print_probe_vector(probes: &Vec<Probe>, probenames_only: bool) {
    if probenames_only {
        for probe in probes {
            println!("{}", probe.name)
        }
    } else {
        println!("name\tstandard_deviation");
        for probe in probes {
            println!("{}", probe)
        }
    }
}

// fn probe_names(probes: &Vec<Probe>) -> Vec<&str> {
//     probes.iter().map(|x| x.name.as_str()).collect()
// }

fn gzip_csv_to_reader(path: &String) -> io::BufReader<flate2::read::GzDecoder<std::fs::File>> {
    // Open File
    let file_result = std::fs::File::open(path);
    let file = match file_result {
        Ok(val) => val,
        Err(err) => panic!("Failed to open file {:?}. Error: {:?}", path, err),
    };

    // Create a reader for the gzipped file
    let decoder = flate2::read::GzDecoder::new(file);
    let gzip_reader: io::BufReader<flate2::read::GzDecoder<std::fs::File>> =
        std::io::BufReader::new(decoder);

    // Return Gzip Reader
    gzip_reader
}

fn get_most_variable_probes(
    mut reader: csv::Reader<std::boxed::Box<dyn std::io::Read>>,
    n_probes_to_return: usize,
    max_proportion_missing_samples: f64,
) -> TopProbes {
    let mut n_probes_failed_qc: usize = 0;
    let mut top_probes: Vec<Probe> = Vec::new();
    let mut min_probe_sd: f64 = 0.0;
    let mut n_total_probes = 0;

    let n_samples: usize = reader.headers().expect("Failed to read headers").len() - 1;

    // For each row in the CSV
    for (rownumber, result) in reader.records().enumerate() {
        n_total_probes += 1;
        let record = match result {
            Ok(val) => val,
            Err(err) => panic!("Failed to read record: {err}"),
        };

        // Get Probe Name
        let probename_opt = record.get(0);
        let probename = match probename_opt {
            Some(val) => val,
            None => panic!("Empty probename at row {}", rownumber),
        };

        // Compute Mean Beta Value for the probe
        let nsamples: usize = record.len() - 1;

        // Get a vector of beta values (skip first element of record since thats the probe name)
        let beta_values: Vec<f64> = record
            .iter()
            .skip(1)
            .filter_map(|x| x.parse::<f64>().ok())
            .collect();

        let n_beta_values = beta_values.len() as f64;
        let missing_beta_values: usize = nsamples - n_beta_values as usize;
        let beta_total: f64 = beta_values.iter().sum();
        let beta_mean = beta_total / n_beta_values;
        let variance: f64 = beta_values
            .iter()
            .map(|b| (b - beta_mean).powf(2.0))
            .sum::<f64>()
            / (n_beta_values - 1.0);
        let sd = variance.sqrt();
        let prop_missing_samples = missing_beta_values as f64 / nsamples as f64 * 100.0;

        // Minumum Probe requirements
        let probe_passes_qc = prop_missing_samples < max_proportion_missing_samples;

        if !probe_passes_qc {
            n_probes_failed_qc += 1
        }

        if top_probes.len() < n_probes_to_return && probe_passes_qc {
            top_probes.push(Probe {
                name: probename.to_string(),
                sd,
            });

            // Sort probes in ascending order (by SD)
            top_probes.sort_by(|x, y| x.sd.partial_cmp(&y.sd).expect("Failed to find min"));
            min_probe_sd = top_probes.first().expect("Lowest SD probe").sd
        } else if sd > min_probe_sd && probe_passes_qc {
            top_probes[0] = Probe {
                name: probename.to_string(),
                sd,
            };

            // Sort probes in ascending order (by SD)
            top_probes.sort_by(|x, y| x.sd.partial_cmp(&y.sd).expect("Failed to find min"));
            min_probe_sd = top_probes.first().expect("Lowest SD probe").sd
        }
    }

    TopProbes {
        top_probes,
        n_probes_failed_qc,
        n_total_probes,
        n_samples,
    }
}

fn eprint_hr() {
    eprintln!("{}", "-".repeat(40));
}

fn get_reader(source: InputSource, delim: u8) -> csv::Reader<std::boxed::Box<dyn std::io::Read>> {
    let rdr = match source {
        InputSource::File(path) => Box::new(gzip_csv_to_reader(&path)),
        InputSource::Stdin => Box::new(io::stdin().lock()) as Box<dyn Read>,
    };

    csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delim)
        .from_reader(rdr)
}

fn stdin_supplied() -> bool {
    !atty::is(atty::Stream::Stdin)
}

fn get_reader_from_stdin_or_file(filepath:Option<&String>, delim:u8) -> csv::Reader<std::boxed::Box<dyn std::io::Read>>{
    let stdin_supplied = stdin_supplied();
    let filepath_supplied = filepath.is_some();
    
    if stdin_supplied && filepath_supplied {
        panic!("Data was supplied to both STDIN and methdata positional argument. Please choose one or the other.")
    };

    let inputsource: InputSource;
    if stdin_supplied {
        inputsource = InputSource::Stdin;
    } else if filepath_supplied {
        let path = filepath.expect("Failed to parse methdata argument");

        inputsource = InputSource::File(path.to_string())
    } else {
        panic!("Must supply methdata either through stdin or methdata positional argument. See `probepicker --help` for details")
    };
    
    eprint_hr();
    eprintln!("Reading Data from {:#?}", inputsource);

    // Return reader
    get_reader(inputsource, delim)
}

fn get_selected_samples(path:String) -> Vec<String>{

    let file = std::fs::File::open(path.as_str()).expect("no such file");
    let buf = std::io::BufReader::new(file);

    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}

fn print_selected_samples(mut reader:csv::Reader<Box<dyn Read>>, samples:Vec<String>){
    
    let headers = reader.headers().expect("Failed to get headers");

    let sample_indices:Vec<usize> = headers
    .iter()
    .enumerate()
    .filter_map(|(i, h)| if samples.contains(&h.to_string()) { Some(i) } else { None })
    .collect();

    let colnames_in_output:Vec<&str> = headers.iter().filter(|x| samples.contains(&x.to_string())).collect();

    let n_samples_found = sample_indices.len();
    let n_samples_searched_for = samples.len();

    eprintln!("{}/{} ({}%) of samples found in file",n_samples_found, n_samples_searched_for, n_samples_found as f64 / n_samples_searched_for as f64 * 100.0);
    
    // Print header
    println!("{}\t{}", headers.get(0).expect("first column name"),colnames_in_output.join("\t"));
    for result in reader.records(){
        let record = result.expect("Failed to parse result");
        let filtered_record: Vec<_> = sample_indices.iter().map(|&i| &record[i]).collect();
        println!("{}\t{}", record.get(0).expect("first column"), filtered_record.join("\t"));
    }
}