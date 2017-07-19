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
	Class to implement alignment of two strings using the Needleman–Wunsch algorithm.
	Runtime is quadratic as we do constant work for each (i,j) for i in s1 and j in s2.

	Class has two public methods: align() and getAlignmentScore()
	"""
	def __init__(self, s1, s2, gapPenalty = GAP_PENALTY, mismatchPenalty = MISMATCH_PENALTY):
		assert type(s1) == type(s2) == type("") # input type must be string
		assert(len(s1)) > 0 # cannot align to an empty string
		assert(len(s2)) > 0
		self._s1 = s1
		self._s2 = s2
		self._gapPenalty = gapPenalty
		self._mismatchPenalty = mismatchPenalty
		self._scoreMatrix = None
		self._bestAlignment = None

	def align(self):
		"""
		Implements Needleman-Wunsch alignment algorithm and returns optimal alignment.

		If align() has previously been called, previous return is memoized and returned
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

	def getScoreMatrix(self):
		return self._scoreMatrix

	def _initializeMatrix(self):
		"""
		Initializes the score matrix for the alignment of s1 and s2

		s1 is placed on the rows of the matrix, while s2 is on the columns

		Helper function for computeScoreMatrix()
		"""
		self._scoreMatrix = [[0 for j in range(len(s2)) + 1] for i in range(len(s1) + 1)]
		self._scoreMatrix[0] = [self._gapPenalty * (j - 1) for j in range(1, len(s2) + 1)]
		for i in range(1, len(s1) + 1):
			self._scoreMatrix[i][0] = self._gapPenalty * (i - 1)


	def _computeScoreMatrix(self):
		"""
		Fills in the score matrix for the alignment of s1 and s2

		Helper function for constructAlignment()
		"""
		self._initializeMatrix()


	def _constructAlignment(self):
		"""
		Uses the completed score matrix
		"""
		if self._scoreMatrix is None:
			self._computeScoreMatrix()


def align(s1, s2):
	"""
	Public function to return an optimal alignment between two nonempty strings s1 and s2.

	Thin wrapper around StringAligner object.
	"""
	aligner = StringAligner(s1, s2)
	return aligner.align()

if __name__ == "__main__":
	pass