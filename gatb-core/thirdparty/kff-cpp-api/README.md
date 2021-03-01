# C++ lib for i/o kff files

Welcome in the C++ library description for Reading or Writing kff files.
This code is not guaranty to be 100% bug-free, so please submit issues if encounter one.

The library is divided into a low level and a high level APIs.
High level API is designed to be easy to use by hiding the section aspect of the file.
The low level API let you read and write sections as you want.
You are entirely free of your kmer organization inside of the file.

The lib is not yet optimized, so it can be slow for big files.
Also, the low level API does not guaranty the file integrity.
It means that if you use the functions in the wrong order, you can write wrong kff files.


# High level reader API

## Enumerating kmers

The high level reader API is made to be very easy to use.
This part of the API allow you to enumerate each pair of kmer / data through the whole file, hiding all the kff datastructures.

The HL API contains 2 main functions:
- *has_next* to know if the enumeration is over or not.
- *next_kmer* that return the next available kmer and its data.
NB: the kmer and data pointer are valid until you use again one of the two functions.

```C++
  Kff_reader reader("test.kff");

  uint8_t * kmer;
  uint8_t * data;

  while (reader.has_next()) {
    reader.next_kmer(kmer, data);
    //[...] use the kmer and its data
  }
```

## Enumerating blocks

If you prefer to enumerate kmers by block instead of one by one, you can use the function *next_block*.

```C++
  Kff_reader reader("test.kff");

  uint8_t * kmers;
  uint8_t * data;

  while (reader.has_next()) {
    uint64_t nb_kmers = reader.next_block(kmers, data);
    //[...] use the kmers and their associated data
  }
```


## How to know properties of my kmers ?

To use the kmers and their data, you need to know some values like k or data_size.
These values of any variable are accessible through the *get_var* method.
You also need to know the encoding used to translate you data from 2-bits to strings.
The *get_encoding* function return an array of the 4 encoded values in order A, C, G, T.

```C++
  // Get the file encoding
  uint8_t * encoding = reader.get_encoding();
  // Get the current k and data_size values
  uint64_t k = reader.get_var("k");
  uint64_t ds = reader.get_var("data_size");
```


# Low level API

## File reader

### Open/close a file

Kff_file is the main object in the library.
This is the object needed to manipulate a binary kff file.

```C++
  // Open a kff file to read it
  Kff_file infile("path/to/file/myfile.kff", "r");
  // [...]
  infile.close();
```

### Header

When you open a file to read it, the software version is automatically read and compared to the file version.
The software version must be equal or greater than the file version.

The encoding of the nucleotides is also automatically read.
It is accessible via the file property *encoding*.

```C++
  // Nucleotide encoding
  uint8_t encoding[4] = infile.encoding;
```

Finally the metadata size is automatically read but not the metadata content.
The size is available via the the file property *metadata_size*.
You can either read the metadata content as follow or let the first call to *read_section_type* skip it for you.

```C++
  // Get the metadate size
  uint32_t size = infile.metadata_size;
  // Allocate memory and get the metadata content
  uint8_t * metadata = new uint8_t[size];
  infile.read_metadata(metadata);
  // Use metadata [...]
  delete[] metadata;
```

### Detect section

Each section start with a special char.
You can use *read_section_type* to get this char and then use the dedicated function to open the corresponding section.

```C++
  char type = infile.read_section_type();
```

### Global variable section

When 'v' char is detected, you can call the the opening of a Section_GV.
The creation of the section on a file will automatically read the variables inside of the section.
The variables are accessible as a std::map<string, uint64_t> with the public section property vars or inside the file property global_vars.

```C++
  // Open the section (automatically read the variables)
  Section_GV sgv(&infile);
  // Read the variables from the section map
  for (auto it : sgv.vars)
    std::cout << it.first << ": " << it.second << std::endl;
  // Read the variables from the global map
  for (auto it : infile.global_vars)
    std::cout << it.first << ": " << it.second << std::endl;
```

### Raw sequences section

When 'r' char is detected, you can call the opening of a Section_Raw.
At the section creation, the number of blocks inside of it is automatically read and stored in the property *nb_blocks*.
The number of remaining blocks in the section is also available in the property *remaining_blocks*.

```C++
  // Get the global variables needed
  uint64_t k = infile.global_vars["k"];
  uint64_t max = infile.global_vars["max"];
  uint64_t ds = infile.global_vars["data_size"];

  // Open the raw section
  Section_Raw sr(&infile);
  cout << sr.remaining_blocks << "/" << sr.nb_blocks << endl;

  // Prepare buffers for sequences and data
  uint8_t * sequence_buffer = new uint8_t[(k + max -1) / 4 + 1];
  uint8_t * data_buffer = new uint8_t[max * ds];

  // Read all the blocks of the section
  for (uint64_t i=0 ; i<sr.nb_blocks ; i++) {
    uint64_t nb_kmers = sr.read_compacted_sequence(sequence_buffer, data_buffer);
    // [...] Use sequences and data
    cout << sr.remaining_blocks << "/" << sr.nb_blocks << endl;
  }

  // Close the raw section
  delete[] sequence_buffer;
  delete[] data_buffer;
  sr.close();
```

**NB**: You can close the section even if you did not read all the blocks of it.
It is triggering a function that will skip the end of the section.


**Warning**: You can directly interact with the file pointer under the hood of Kff_file objects.
If you do so, variables like *remaining_blocks* and automatic readings procedures on section closing can become wild.
So, please don't rely on them after a direct file pointer usage or manually update them.


### Minimizer sequences section

When 'm' char is detected, you can call the the opening of a Section_Minimizer.
This section is very similar from the raw section.
In fact, you can use the exact same code after the opening.
But where raw blocks reading function uses direct file reading, the one here hide more computation.
This computation overhead is due to the fact that minimizer are not stored inside of the blocks but at the beginning of the section.

To avoid the computation, you can directly read sequences using minimizer section specialized function that read the sequence without the minimizer.
The minimizer is accessible in the section attribute *minimizer*.

```C++
  // [...] Same previous variables
  uint64_t m = infile.global_vars["m"];

  // Open the minimizer section
  Section_Minimizer sm(infile);
  uint8_t * mini = new uint8_t[m / 4 + 1];
  memcpy(mini, sm.minimizer, m%4==0 ? m/4 : m/4+1);
  // [...] Use the minimizer

  // Read all the sequences
  for (uint64_t i=0 ; i<sm.nb_blocks ; i++) {
    uint64_t minimizer_position;
    uint64_t nb_kmers = sm.read_compacted_sequence_without_mini(sequence_buffer, data_buffer, minimizer_position);
    // [...] Use sequence and data
    cout << "minimizer position: " << minimizer_position << endl;
  }

  // Close the section (and detroy the section minimizer memory!)
  sm.close();
  delete[] mini;
```

**Warning**: The minimizer memory is deleted on section closing.
So, if needed, remember to copy the value before.


### Reading section

Because Raw sections and Minimizer sections share a lot of their design for reading functions, they both inherit from a class called *Block_section_reader*.
This class contains all the properties and functions described in the raw sequence section part.

These sections are created using a static function that returns a nullptr if the file cannot recognize a raw or minimizer section.

```C++
  // Construct a block reader
  Block_section_reader * br = Block_section_reader::construct_section(&infile);

  if (br == nullptr) {
    std::cerr("Next section is not a sequence section");
  } else {
    // [...] Do stuff here
  }
```

### Skipping blocks or sections

For some reasons, you may want to jump over a block inside of a reading section or over a complete section.
For that, we create a dedicated function for each:

```C++
  // Construct a block reader
  Block_section_reader * br = Block_section_reader::construct_section(&infile);
  // Skip a block inside of a reader
  br.jump_sequence();
  br.close();

  // jump over a complete section
  bool jumped = infile.jump_next_section();
  // [...] jumped == false if no section can be jumped
```


## File writer

### Open/close a file

To open a kff file in writing mode, just put the 'w' flag at construction.
NB: All the directories of the path to the file that you want to create must already exist.

```C++
  // Open a kff file to write inside of it.
  // Path/to/file/ must exist prior to opening.
  Kff_file outfile("path/to/file/myfile.kff", "w");
  // [...]
  outfile.close();
```

### Header

The header first contains the kff file version. This 2 Bytes value is automatically added by the API at the beginning of the file.
Then you have to write the encoding that you use using one of the *write_encoding* functions.
If you want, you can finish the header by adding some metadata (Byte array format).
If you do not set any metadata, the API will automatically create a 0 Byte metadata field.

```C++
  //                    A  C  G  T
  uint8_t encoding[] = {0, 1, 3, 2};
  // Encoding writing alternative 1
  outfile.write_encoding(encoding);
  // Encoding writing alternative 2
  // outfile.write_encoding(0, 1, 3, 2);

  // Write metadata
  std::string meta_str = "<3 KFF <3";
  outfile.write_metadata(meta_str.length(), (uint8_t *)meta.c_str());
```


### Global Variable section

GV sections start with the number of variables described in it.
Initially this value is automatically set to 0 and updated when you call the close function.

```C++
  // Open the section
  Section_GV sgv(outfile);
  
  // Write all the values needed
  sgv.write_var("k", 7);
  sgv.write_var("m", 3);
  sgv.write_var("max", 8);
  sgv.write_var("data_size", 1);

  sgv.close();
```


### String sequences to binary

This API does not provide any function to transform a DNA string into a byte array coding for the binary sequence.
To see how to properly encode a sequence, please refer to the kmer file format reference.

For the following parts, we assume that your code include a function *encode* that take a string and return a uint8_t \* with the correct format.


### Raw sequences section

As in the GV section, the raw section starts with a number of blocks that is automatically set by the API on closing.
To write raw sections, you just have to use the *write_compacted_sequence*.
This function take the byte array for your 2-bits nucleotidic sequence and a byte array for their related data.

Prior to raw section writings, do not forget to write the variables needed (cf kff reference).

```C++
  // Open the raw section
  Section_Raw sr(outfile);

  std::string sequence;
  uint8_t counts[8];
  uint8_t * byte_seq;

  // Encode one kmer
  sequence = "GATTACA";
  uint8_t * byte_seq = encode(sequence);
  counts[0] = 3;
  // Write the 7-mer
  sr.write_compacted_sequence(byte_seq, sequence.length(), counts);

  // Encode 8 overlapping 7-mers: TGGTACA, GGTACAA, GTACAAG, ...
  sequence = "TGGTACAAGTTACC";
  counts[0]=1; counts[1]=221; counts[2]=7; counts[3]=1;
  counts[4]=5; counts[5]=6; counts[6]=12; counts[7]=198;
  sr.write_compacted_sequence(byte_seq, sequence.length(), counts);

  sr.close();
```


### Minimizer sequences section

In the minimizer section, everything is very similar to the raw sequences.
At the beginning of the section, you need to write the minimizer using the method *write_minimizer*.
To write a block, you only need one more piece of information than for a row block, which is the index where the minimizer is present in the sequence.

```C++
  // Open the raw section
  Section_Minimizer sm(outfile);

  std::string sequence;
  uint8_t counts[8];
  uint8_t * byte_seq;

  // Write the minimizer
  byte_seq = encode("ACA");
  sm.write_minimizer(byte_seq);

  // Encode one kmer
  sequence = "GATTACA";
  //              ^ minimizer (pos = 4)
  uint8_t * byte_seq = encode(sequence);
  counts[0] = 3;
  // Write the 7-mer
  sm.write_compacted_sequence(byte_seq, sequence.length(), 4, counts);
  //                                    minimizer position ^

  // Encode 8 overlapping 7-mers: TGGTACA, GGTACAA, GTACAAG, ...
  sequence = "GGTACAAGTTACCT";
  //             ^ minimizer (pos = 3)
  counts[0]=1; counts[1]=221; counts[2]=7; counts[3]=1;
  counts[4]=5; counts[5]=6; counts[6]=12; counts[7]=198;
  sm.write_compacted_sequence(byte_seq, sequence.length(), 3, counts);

  sm.close();
```

