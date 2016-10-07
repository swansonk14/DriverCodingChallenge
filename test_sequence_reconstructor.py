import random
import sys
import os
import sequence_reconstructor

DNA_letters = 'ACTG'

# Given the number of segments desired and the length
# of a segment, generate and return a test sequence
# segments for that sequence.
def test_generator(num_segments, segment_length):
	# Initialize variables
	sequence = ''
	segments = []
	segment = ''

	# Generate first segment
	for i in range(segment_length):
		letter = random.choice(DNA_letters)
		segment += letter

	sequence += segment
	segments.append(segment)

	# Generate remaining segments so that each
	# overlaps with the previous by some random
	# amount between 50% and 100%
	for i in range(1, num_segments):
		index = random.randint(0, segment_length / 2 - 1)
		# Overlap between beginning of this segment and
		# end of previous segment
		segment = segments[i-1][index:]

		# Fill in the remaining letters
		for j in range(segment_length - len(segment)):
			letter = random.choice(DNA_letters)
			segment += letter
			sequence += letter

		segments.append(segment)

	return sequence, segments

# Convert the segments to FASTA form
# and save to a file
def segments_to_file(segments, filename):
	with open(filename, 'wb+') as segments_file:
		for i in range(len(segments)):
			segment = segments[i]
			segments_file.write('>FRAG_' + str(i) + '\n')
			segments_file.write(segment + '\n')

# Read in reconstructed sequence from file
def read_sequence(filename):
	with open(filename, 'rb+') as sequence_file:
		return sequence_file.read()

# Run the test
def run_test(num_segments, segment_length):
	segments_file = 'test_sequence_reconstructor_segments.txt'
	sequence_file = 'test_sequence_recontructor_sequence.txt'

	# Generate and save segments
	print "Generating sequence..."
	sequence, segments = test_generator(num_segments, segment_length)
	segments_to_file(segments, segments_file)

	# Run the sequence reconstructor
	print "Reconstructing sequence..."
	sequence_reconstructor.run_reconstructor(segments_file, sequence_file)

	# Get reconstructed sequence
	reconstructed_sequence = read_sequence(sequence_file)

	# Remove temporary files
	os.remove(segments_file)
	os.remove(sequence_file)

	if sequence == reconstructed_sequence:
		print "The sequence reconstructor was successful!"
	else:
		print "The sequence reconstructor failed."
	print "Correct sequence:\n" + sequence
	print "Reconstructed sequence:\n" + reconstructed_sequence


#--------------------------------------------------------


# Main function
if __name__ == '__main__':
	# Read in arguments
	if len(sys.argv) < 3:
		exit("Error parsing input. Usage: \"python test_sequence_reconstructor.py <number_of_segments> <segment_length>\"")
	num_segments = int(sys.argv[1])
	if num_segments <= 1:
		exit("Error, please select more than 1 segment")
	segment_length = int(sys.argv[2])
	if segment_length <= 2:
		exit("Error, please select a segment length greater than 2")

	run_test(num_segments, segment_length)