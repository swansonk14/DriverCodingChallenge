# Coding challenge for Driver Group internship.

## Running the code

To run the code, type:

```
python sequence_reconstructor.py <file_with_segments> <file_for_reconstructed_sequence>
```

where &lt;file_with_segments&gt; is the file containing the DNA segments in FASTA format and &lt;file_for_reconstructed_sequence&gt; is the file in which you wish to save the reconstructed sequence.

For example, to reconstruct the sequence provided in the coding challenge, type:

```
python sequence_reconstructor.py dataset.txt dataset_solution.txt
```

and then view dataset_solution.txt to see the reconstructed sequence.

## Testing the code

To test the code, type:

```
python test_sequence_reconstructor.py <number_of_segments> <segment_length>
```

where &lt;number_of_segments&gt; is the number of DNA segments you wish to use (at least 1) and &lt;segment_length&gt; is the length of a segment (at least 3). This will generate a random DNA sequence and will break it down randomly into segments that overlap by at least 50%. It will then run the reconstructor and check that the reconstruction matches the original sequence.

## Explanation of the solution

Please see README.pdf for an explanation of how the reconstruction is done.
