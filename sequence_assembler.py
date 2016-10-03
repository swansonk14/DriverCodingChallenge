# Read in a file in FASTA format and return a
# dictionary of IDs to DNA segments.
def read_segments(filename):
	with open(filename, 'rb+') as dna_file:
		segment_dict = {}
		segment = ''
		segment_id = ''
		for line in dna_file:
			# If line is the beginning of a new segment,
			# get ID and initialize dictionary element.
			if line[0] == '>':
				segment_id = line[1:].strip()
				segment_dict[segment_id] = ''
			# Combine lines of a DNA segment into one string.
			else:
				segment_dict[segment_id] += line.strip()
	return segment_dict

# Takes in a dictionary of IDs to segments and returns
# a dictionary mapping IDs to an array of tuples, where
# each tuple contains the ID of a segment that overlaps
# with more than 50% of this segment along with the index
# of this segment at which the overlap begins. 
def find_overlaps(segment_dict):
	overlap_dict = {}
	# Initialize overlap dictionary with IDs.
	for segment_id in segment_dict.keys():
		overlap_dict[segment_id] = []
	# A dictionary mapping a subsequence of DNA
	# to an array of tuples of segment ID and index
	# at which the subsequence starts.
	subsequence_to_matches = {}
	# Loop through segments
	for segment_id in segment_dict.keys():
		print segment_id
		segment = segment_dict[segment_id]
		# Check subsequences starting at the beginning of the segment
		# which are at least 50% of the length of the segment.
		for index in range(len(segment) / 2, len(segment)):
			subsequence = segment[:index]
			if not subsequence in subsequence_to_matches.keys():
				subsequence_to_matches[subsequence] = [(segment_id, 0)]
			else:
				# Identify overlaps
				for overlap in subsequence_to_matches[subsequence]:
					overlap_id = overlap[0]
					overlap_start = overlap[1]
					# If the matched segment ends with this subsequence
					# or if the two segments are identical, then add this
					# segment as a successor to the matched segment.
					if (overlap_start > 0 or (overlap_start == 0 and index == len(segment) - 1)):
						overlap_dict[overlap_id].append((segment_id, overlap_start))
				subsequence_to_matches[subsequence].append((segment_id, 0))
		# Check subsequences ending at the end of the segment
		# which are at least 50% of the length of the segment.
		for index in range(len(segment) / 2 + 1):
			subsequence = segment[index:]
			if not subsequence in subsequence_to_matches.keys():
				subsequence_to_matches[subsequence] = [(segment_id, index)]
			else:
				# Identify overlaps
				for overlap in subsequence_to_matches[subsequence]:
					overlap_id = overlap[0]
					overlap_start = overlap[1]
					# If the matched segment starts with this subsequence,
					# then add the matched segment as a successor to this segment.
					if (overlap_start == 0):
						overlap_dict[segment_id].append((overlap_id, index))
				subsequence_to_matches[subsequence].append((segment_id, index))
	return overlap_dict

# Take dictionaries from IDs to segments and from IDs
# to overlaps and combine the segments into the
# original sequence.
def combine_segments(segment_dict, overlap_dict):
	sequence = ''
	# stack = []
	# num_segments = len(segment_dict.keys())
	# segments_seen = set()
	# # Try each segment as starting segment
	# for segment in segment_dict.keys():
	# 	stack.append((segment_id))
	# 	segments_seen.add(segment_id)


	return sequence

# Save the sequence to the given file.
def save_sequence(sequence, filename):
	with open(filename, 'wb+') as sequence_file:
		sequence_file.write(sequence)


# Main
segments_file = 'dataset.txt'
segment_dict = read_segments(segments_file)
# overlap_dict = find_overlaps(segment_dict)
# sequence = combine_segments(segment_dict, overlap_dict)
# sequence_file = 'sequence.txt'
# save_sequence(sequence, sequence_file)