"""
Package to find the optimal alignment between two strings s1 and s2.

Runtime of alignment is quadratic in the length of the input strings.
"""

# We need the penalty to be greater than the penalty for a gap, but strictly less than double
# the penalty for a gap. If the mismatch penalty were <= 2 * the gap penalty, we would get
# unreasonable alignments (e.g. for s1 = AGT and s2 = ACT we might get A-GT and AC-T)

# For biological sequences, an even better scoring system would take into account that deletions
# tend to have more than one base, so the first gap should be penalized more heavily than subsequent
# gaps. However, this is okay for a fist approximation.
GAP_PENALTY = 2
MISMATCH_PENALTY = 3

class StringAligner(object):
	"""
	Class to implement alignment of two strings using the Needleman-Wunsch algorithm.
	Runtime is quadratic as we do constant work for each (i,j) for i in s1 and j in s2.

	N.B. StringAligner is NOT case-sensitive! "A" is not the same as "a" and will be considered
	a mismatch.

	gapPenalty and mismatchPenalty should be greater than 0, since this implementation of the
	algorithm seeks to minimize the alignment score.

	Class has two public methods: align() and getAlignmentScore()
	It also implements a __repr__ method so the object can be printed
	"""
	def __init__(self, s1, s2, gapPenalty = GAP_PENALTY, mismatchPenalty = MISMATCH_PENALTY):
		assert type(s1) == type(s2) == type("") # input type must be string
		assert(gapPenalty > 0)
		assert(mismatchPenalty > 0)
		self._s1 = s1
		self._s2 = s2
		self._gapPenalty = gapPenalty
		self._mismatchPenalty = mismatchPenalty
		self._scoreMatrix = None
		self._bestAlignment = None

	def __repr__(self):
		"""
		Allow for a nice print representation of the optimal alignment.
		"""
		if self._bestAlignment is not None:
			return "s1: {0} \n s2: {1}".format(self._bestAlignment[0], self._bestAlignment[1])
		return ""

	def align(self):
		"""
		Implements Needleman-Wunsch alignment algorithm and returns optimal alignment.

		If align() has previously been called, previous return is memoized and returned

		Return type is a tuple of (alignment for s1, alignment for s2)
		"""
		if self._bestAlignment is None:
			self._constructAlignment()
		return self._bestAlignment

	def getAlignmentScore(self):
		"""
		Returns alignment score of best alignment found for s1 and s2
		"""
		if self._scoreMatrix is None:
			self._computeScoreMatrix()
		return self._scoreMatrix[-1][-1]

	def _initializeMatrix(self):
		"""
		Initializes the score matrix for the alignment of s1 and s2

		s1 is placed on the rows of the matrix, while s2 is on the columns

		Helper function for computeScoreMatrix()
		"""
		# initialize an len(s1) + 1 x len(s2) + 1 matrix of 0s
		self._scoreMatrix = [[0 for j in range(len(self._s2) + 1)] for i in range(len(self._s1) + 1)]

		# set the first row and first column to index * gapPenalty
		self._scoreMatrix[0] = [self._gapPenalty * j for j in range(len(self._s2) + 1)]
		for i in range(len(self._s1) + 1):
			self._scoreMatrix[i][0] = self._gapPenalty * i


	def _computeScoreMatrix(self):
		"""
		Fills in the score matrix for the alignment of s1 and s2

		Helper function for constructAlignment()
		"""
		self._initializeMatrix()
		for i in range(1, len(self._s1) + 1):
			for j in range(1, len(self._s2) + 1):
				
				# mismatch score for this part of the alignment is 0 if the characters match in s1 and s2; otherwise it is mismatchPenalty
				mismatch = int(self._s1[i - 1] != self._s2[j - 1]) * self._mismatchPenalty # matrix indices are off by one from string indices, so we have to do i - 1 and j - 1

				# recurrence relation: we have three options for obtaining the minimum score at (i,j)
				# 1) traverse only s1 (i.e. take the score at (i - 1, j) and add the gap penalty
				# 2) traverse only s2 (i.e. take the score at (i, j - 1) and add the gap penalty)
				# 3) traverse both s1 and s2 (i.e. take the score at (i - 1, j - 1) and add a mismatch penalty if the characters don't match)
				self._scoreMatrix[i][j] = min(
					self._scoreMatrix[i - 1][j] + self._gapPenalty
					, self._scoreMatrix[i][j - 1] + self._gapPenalty
					, self._scoreMatrix[i - 1][j - 1] + mismatch
				)


	def _constructAlignment(self):
		"""
		Uses the completed score matrix to return the optimal alignment between s1 and s2

		Alignments are stored as a tuple of (alignment for s1, alignment for s2)
		"""
		if self._scoreMatrix is None:
			self._computeScoreMatrix()

		i, j = len(self._s1), len(self._s2) # keep track of s1 and s2
		s1Align, s2Align = "", ""

		while i > 0 or j > 0:
			if i == 0:
				# we are at the beginning of s1
				s2Align += self._s2[j - 1]
				s1Align += "-"
				j -= 1
			elif j == 0:
				# we are at the beginning of s2
				s1Align += self._s1[i - 1]
				s2Align += "-"
				i -= 1
			else: 
				# we still need to traverse both s1 and s2
				localOptimum = min(self._scoreMatrix[i - 1][j], self._scoreMatrix[i][j - 1], self._scoreMatrix[i - 1][j - 1])
				if localOptimum == self._scoreMatrix[i - 1][j]:
					s1Align += self._s1[i - 1]
					s2Align += "-"
					i -= 1
				elif localOptimum == self._scoreMatrix[i][j - 1]:
					s2Align += self._s2[j - 1]
					s1Align += "-"
					j -= 1
				else:
					s1Align += self._s1[i - 1]
					s2Align += self._s2[j - 1]
					i -= 1
					j -= 1

		# since we moved from the end of the matrix to the source, our alignments are backward
		s1Align = s1Align[::-1] # hacky slice syntax for reversing strings
		s2Align = s2Align[::-1]

		self._bestAlignment = (s1Align, s2Align)

def align(s1, s2):
	"""
	Public function to return an optimal alignment between two nonempty strings s1 and s2.

	Thin wrapper around StringAligner object.
	"""
	aligner = StringAligner(s1, s2)
	return aligner.align()

if __name__ == "__main__":
	# example code
	S1 = "AGGCT"
	S2 = "AGCA"
	aligner = StringAligner(S1, S2)
	aligner.align()
	print aligner
	print align("xyxyx", "")
