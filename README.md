# DNA Bloom Filter
This is a complete code flow of Bloom Filter application in DNA storage. It can achieve two functions, one is the 
anti-contamination function based on the Bloom Filter (BF), another is the file version control function based on the 
Counting Bloom Filter (CBF). 
It includes the generation and detection of the BF and CBF, the simulation of the synthesizing, sequencing and error 
sequences generation. 

## Environment Configuration
The kit is developed by **Python3.9**.

In addition, the packages we are calling now is as follows:

- [x] numpy
- [x] pandas
- [x] pickle
- [x] collections
- [x] sys
- [x] os
- [x] random
- [x] math
- [x] copy
- [x] datetime
- [x] unittest
- [x] matplotlib


## Kit Tree Diagram
```html
├── bean                               // BF and CBF module
│    ├── bloom.py                      // Training and checking process of the BF and CBF
│    ├── hasher.py                     // Hash functions calculations
├── examples                            // Examples of running scripts
│    ├── simulation_actual_data.py     // The simulation process and anti-contamination function for the actual data
│    ├── simulation_random_data.py     // The simulation process and anti-contamination function for the random data
│    ├── simulation_file_version.py    // The simulation process and file file version control function for different versions files
├── exps                               // The process of the experiment
│    ├── operations.py                 // The simulation of synthesizing, sequencing and error sequences generation
│    ├── pipeline.py                   // Main calling function
├── files                              // The directory of the input files and output files
│    ├── original files                // The original files encoded as DNA sequences used for experimental analysis
│    ├── version                       // Different versions of the files and their corresponding encoded DNA files
│    ├── output                        // The directory of the output files
├── tests                              // Test module
│    ├── test_delete.py                // Test for BF detection of the sequences with deletion errors
│    ├── test_insert.py                // Test for BF detection of the sequences with insertion errors
│    ├── test_replace.py               // Test for BF detection of the sequences with replacement errors
├── utils                              // Utils module
│    ├── log.py                        // Logs outputting in console
│    ├── model_saver.py                // Save the model to the .pkl file, or load the model from the .pkl file
│    ├── recorder.py                   // Record the variable during synthesizing
│    ├── recorder_seq.py               // Record the variable in a certain depth during sequencing
├── README.md                          // Description document
```

## Parameters
- **expected_rate[float]** The target FPR of the BF or CBF
- **payload_length[int]** The payload length of the DNA sequences
- **copy_number[int]** The molecular copy number of each DNA sequence
- **stdev[int]** The standard deviation of normal distribution after synthesis of DNA sequences
- **hash_function_type[str]** The hash type for the generation of the BF or CBF
- **read_depth[int]** The read depth of the sequencing
- **root_path[str]** The root path of the files
- **insertion_rate[float]** The insertion rate for the error sequences generation, default value is 0.00075
- **deletion_rate[float]** The deletion rate for the error sequences generation, default value is 0.00075
- **substitution_rate[float]** The substitution rate for the error sequences generation, default value is 0.0015
- **seed[int]** The random seed, default value is 30
- **actual_dnas[bool]** Whether to use the actual DNA sequences
- **digital_count[int]** The number of DNA sequences need to be randomly generated
- **training_dnas[list]** Actual DNA sequences for training the BF or CBF
- **checking_dnas[list]** Actual DNA sequences for checking the BF or CBF
- **training_name[str]** The name of the training_dnas
- **checking_name[str]** The name of the checking_dnas
- **filter_size[int]** The expected size of the array from BF or CBF. The defalut value is 0, which means it will be 
                       calculated based on the target FPR of the BF or CBF and the number of the DNA sequences
- **hash_size[int]** The number of hash functions for the BF and CBF generation
- **counting_type[bool]** Whether need to generate a CBF
- **record_coverage[bool]** Whether to record the number of the checked DNA sequences checked in each coverage
- **coverage_threshold[int]** The coverage threshold for removing false positive sequences
- **need_record[bool]** Whether record sequencing data and purely random sequences
- **verbose[bool]** Whether need to show process

## Running
Please see the running scripts in examples folder.
