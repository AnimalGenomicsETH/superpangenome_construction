// faster path concordance
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        println!("Usage script file1 file2");
        process::exit(1);
    }
    let mut args = args.iter();

    args.next();

    let f1 = fs::File::open(args.next().unwrap()).expect("File not found!");
    let f2 = fs::File::open(args.next().unwrap()).expect("File not found!");

    let mut reader1 = BufReader::new(f1).lines().map(|x| x.unwrap());
    let mut reader2 = BufReader::new(f2).lines().map(|x| x.unwrap());

    fn get_vector(line: &str) -> Vec<u32> {
        line.to_string()
            .split(",")
            .map(|x| x.parse::<u32>().unwrap())
            .collect()
    }

    // ID of the animal

    let mut int1: u32;
    let mut chrom1: u32;
    let mut vec1: Vec<u32>;
    let mut int2: u32;
    let mut chrom2: u32;
    let mut vec2: Vec<u32>;
    let mut anim1_list: Vec<String> = Vec::new();
    let mut anim2_list: Vec<String> = Vec::new();
    let mut total_geno: u32 = 0;
    let mut concor_geno: u32 = 0;

    loop {
        let header: String = reader1.next().unwrap();

        if header.starts_with("#CHROM") {
            anim1_list = header.split(",").map(|x| x.to_string()).collect();
            anim1_list = anim1_list[2..].to_vec();
            //println! {"{:?}",anim1_list};
        } else {
            vec1 = get_vector(&header);
            chrom1 = vec1[0];
            int1 = vec1[1];
            break;
        }
    }

    loop {
        let header: String = reader2.next().unwrap();

        if header.starts_with("#CHROM") {
            anim2_list = header.split(",").map(|x| x.to_string()).collect();
            anim2_list = anim2_list[2..].to_vec();
        } else {
            vec2 = get_vector(&header);
            chrom2 = vec2[0];
            int2 = vec2[1];
            break;
        }
    }

    println! {"{:?}",anim1_list};

    loop {
        if chrom1 == chrom2 {
            if int1 == int2 {
                // let dict_geno1: HashMap<_, _> = anim1_list.iter().zip(anim2_list.iter()).collect();
                let dict_geno1: HashMap<_, _> = anim1_list.iter().zip(vec1[2..].iter()).collect();
                let dict_geno2: HashMap<_, _> = anim2_list.iter().zip(vec2[2..].iter()).collect();
                println!(
                    "Overlap at {} {} {} {:?} {:?} {:?}",
                    chrom1, chrom2, int1, anim1_list, dict_geno1, dict_geno2
                );
                for (key, value) in dict_geno1.iter() {
                    // println!("key {} value {}", key, value);
                    // println!("{} {:?}", key, dict_geno2.get(key));
                    match dict_geno2.get(key) {
                        Some(x) => {
                            total_geno += 1;
                            if x == value {
                                concor_geno += 1
                            }
                        }
                        _ => {}
                    }
                }
                vec1 = match reader1.next() {
                    Some(x) => get_vector(&x),
                    None => break,
                };
                chrom1 = vec1[0];
                int1 = vec1[1];
            }
            if int1 < int2 {
                vec1 = match reader1.next() {
                    Some(x) => get_vector(&x),
                    None => break,
                };
                chrom1 = vec1[0];
                int1 = vec1[1];
            } else if int1 > int2 {
                vec2 = match reader2.next() {
                    Some(x) => get_vector(&x),
                    None => break,
                };
                chrom2 = vec2[0];
                int2 = vec2[1];
            }
        } else {
            println!("{} {}", chrom1, chrom2);
            if chrom1 < chrom2 {
                loop {
                    vec1 = match reader1.next() {
                        Some(x) => get_vector(&x),
                        None => break,
                    };
                    chrom1 = vec1[0];
                    int1 = vec1[1];
                    if chrom1 == chrom2 {
                        break;
                    }
                }
            }

            if chrom1 > chrom2 {
                loop {
                    vec2 = match reader2.next() {
                        Some(x) => get_vector(&x),
                        None => break,
                    };
                    chrom2 = vec2[0];
                    int2 = vec2[1];
                    if chrom1 == chrom2 {
                        break;
                    }
                }
            }
        }
    }

    let concor_rate: f32 = concor_geno as f32 * 100.00 / total_geno as f32;

    println!("{} {} {:.2}", total_geno, concor_geno, concor_rate);
}
