import random

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
		letter = choice(DNA_letters)
		segment += letter

	sequence += segment
	segments.append(segment)

	# Generate remaining segments so that each
	# overlaps with the previous by some random
	# amount between 50% and 100%
	for i in range(num_segments):
		index = random.randint(0, segment_length / 2 - 1)
		# Overlap between beginning of this segment and
		# end of previous segment
		segment = segments[i][index:]

		# Fill in the remaining letters
		for j in range(segment_length - len(segment)):
			segment += choice(DNA_letters)

		sequence += segment
		segments.append(segment)

	return sequence, segments


#--------------------------------------------------------


# Main function
if __name__ == '__main__':
	# Read in arguments
	if len(sys.argv) < 3:
		exit("Error parsing input. Usage: \"python sequence_assembler.py <file_with_segments> <file_for_reconstructed_sequence>\"")
	segments_file = sys.argv[1]
	sequence_file = sys.argv[2]

	# Perform sequence reconstruction
	segment_dict = read_segments(segments_file)
	overlap_dict = find_overlaps(segment_dict)
	print overlap_dict
	sequence = combine_segments(segment_dict, overlap_dict)
	save_sequence(sequence, sequence_file)

	print "Done!"