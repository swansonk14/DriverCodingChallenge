import sys

# Read in a file in FASTA format and return a dictionary
# with DNA segment IDs mapping to DNA segments
def read_segments(filename):
	with open(filename, 'rb+') as dna_file:
		segment_dict = {}
		segment_id = ''
		segment = ''
		for line in dna_file:
			# If a line starts with '>', it is the
			# beginning of a new segment
			if line[0] == '>':
				segment_id = line[1:].strip()
				segment_dict[segment_id] = ''
			# Otherwise it is part of the current segment
			else:
				segment_dict[segment_id] += line.strip()
	return segment_dict

# Takes in a dictionary mapping segment IDs to segments.
# Returns a dictionary mapping IDs to a dictionary of
# overlapping segments. This overlapping dictionary maps
# the overlapping segment's ID to the index in the current
# segment at which the overlap begins. Only overlaps of over
# 50% of the segment are considered
# Ex. If segment a is ACT, segment b is CTG, segment c
# is TAAG, and segment d is ACTGG then we get: {a: {b: 1, d:0}
def find_overlaps(segment_dict):
	# Initialize overlap dictionary with all segment IDs
	overlap_dict = {}
	for segment_id in segment_dict.keys():
		overlap_dict[segment_id] = {}

	# Initialize a dictionary which maps a subsequence of DNA
	# to an array of tuples of the form (segment_id, index), where
	# segment_id is the ID of a segment which contains this subsequence
	# and index is the index in the segment at which the subsequence begins.
	# This allows for comparisons between subsequences of all the segments
	# simultaneously.
	subsequence_to_matches = {}

	# Loop through segments and find overlaps
	for segment_id in segment_dict.keys():
		segment = segment_dict[segment_id]

		# Check for overlaps with subsequences of over 50% of the length
		# of the segment starting at the beginning of the segment
		for index in range(len(segment) / 2 + 1, len(segment)):
			subsequence = segment[:index]

			# If the subsequence has not been seen before,
			# add it to our dictionary of subsequence
			if not subsequence in subsequence_to_matches:
				subsequence_to_matches[subsequence] = [(segment_id, 0)]

			# If the subsequence has been seen before,
			# then identify all the other segments with
			# the same subsequence. These are the segments
			# which overlap with our current segment.
			else:
				for overlap in subsequence_to_matches[subsequence]:
					overlap_id = overlap[0]
					overlap_start = overlap[1]

					# If the matched segment ends with this subsequence
					# or if the two segments are identical, then add this
					# segment as a successor to the matched segment.
					if (overlap_start > 0 or (overlap_start == 0 and index == len(segment) - 1)):
						overlap_dict[overlap_id][segment_id] = overlap_start
				subsequence_to_matches[subsequence].append((segment_id, 0))

		# Check for overlaps with subsequences of at least 50% of the length
		# of the segment ending at the end of the segment
		for index in range(1, len(segment) / 2 + 1):
			subsequence = segment[index:]

			# If the subsequence has not been seen before,
			# add it to our dictionary of subsequence
			if not subsequence in subsequence_to_matches:
				subsequence_to_matches[subsequence] = [(segment_id, index)]

			# If the subsequence has been seen before,
			# then identify all the other segments with
			# the same subsequence. These are the segments
			# which overlap with our current segment.
			else:
				for overlap in subsequence_to_matches[subsequence]:
					overlap_id = overlap[0]
					overlap_start = overlap[1]
					# If the matched segment starts with this subsequence,
					# then add the matched segment as a successor to this segment.
					if (overlap_start == 0):
						overlap_dict[segment_id][overlap_id] = index
				subsequence_to_matches[subsequence].append((segment_id, index))

	return overlap_dict

# Takes dictionaries from IDs to segments and from IDs to overlaps.
# Reconstructs and returns the original DNA sequence.
def combine_segments(segment_dict, overlap_dict):
	# Try each segment as the starting segment of the sequence
	for starting_segment_id in segment_dict.keys():
		# Initailize variables
		stack = []
		depth = 0
		previous_depth = -1
		previous_segment_id = None
		num_segments = len(segment_dict.keys())
		segments_seen = set()
		segment_order = []

		# Initialize this segment as first node
		node = (starting_segment_id, depth)
		stack.append(node)

		# Build a sequence starting with the above node
		while len(stack) > 0 and depth != num_segments - 1:
			# Get next segment
			node = stack.pop()
			current_segment_id = node[0]
			depth = node[1]

			# If we are backtracking, delete previous segment_id
			# since it's no longer in the sequence we're constructing
			while depth <= previous_depth:
				previous_segment_id = segment_order.pop()
				segments_seen.remove(previous_segment_id)
				previous_depth -= 1

			# Add current node to segments seen and segment order
			segments_seen.add(current_segment_id)
			segment_order.append(current_segment_id)
			previous_depth = depth

			# Add successor segments
			for next_segment_id in overlap_dict[current_segment_id].keys():
				# Skip segments already in our constructed sequence
				if next_segment_id in segments_seen:
					continue
				# Add successor
				node = (next_segment_id, depth + 1)
				stack.append(node)

		# Check if we found the sequence
		if depth == num_segments - 1:
			break

	# Construct sequence by looping through the segments
	# in order and adding each segment up to the point it
	# overlaps with the next segment
	sequence = ''
	for i in range(len(segment_order) - 1):
		segment_id = segment_order[i]
		next_segment_id = segment_order[i + 1]
		overlap_index = overlap_dict[segment_id][next_segment_id]
		segment = segment_dict[segment_id]
		sequence += segment[:overlap_index]

	# Add final segment
	segment_id = segment_order[len(segment_order) - 1]
	segment = segment_dict[segment_id]
	sequence += segment

	return sequence

# Save the sequence to the given file
def save_sequence(sequence, filename):
	with open(filename, 'wb+') as sequence_file:
		sequence_file.write(sequence)

# Run the sequence reconstructor
def run_reconstructor(segments_file, sequence_file):
	# Perform sequence reconstruction
	segment_dict = read_segments(segments_file)
	overlap_dict = find_overlaps(segment_dict)
	sequence = combine_segments(segment_dict, overlap_dict)
	save_sequence(sequence, sequence_file)


#--------------------------------------------------------


# Main function
if __name__ == '__main__':
	# Read in arguments
	if len(sys.argv) < 3:
		exit("Error parsing input. Usage: \"python sequence_assembler.py <file_with_segments> <file_for_reconstructed_sequence>\"")
	segments_file = sys.argv[1]
	sequence_file = sys.argv[2]

	run_reconstructor(segments_file, sequence_file)

	print "Done!"