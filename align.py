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
GAP_OPENING_PENALTY = 2
GAP_EXTENSION_PENALTY = 1
MISMATCH_PENALTY = 1
MATCH_BONUS = -2

def zeroMatrix(nRow, nCol):
	return [[0 for j in range(nCol)] for i in range(nRow)]

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
	def __init__(self, s1, s2, gapOpeningPenalty = GAP_OPENING_PENALTY, gapExtensionPenalty = GAP_EXTENSION_PENALTY, mismatchPenalty = MISMATCH_PENALTY, matchBonus = MATCH_BONUS):
		assert type(s1) == type(s2) == type("") # input type must be string
		assert gapOpeningPenalty > 0
		assert gapExtensionPenalty >= 0 # if 0 we use linear scoring; if greater than 0 we use affine scoring
		assert mismatchPenalty > 0
		assert matchBonus <= 0
		self._s1 = s1
		self._s2 = s2
		self._gapOpeningPenalty = gapOpeningPenalty
		self._gapExtensionPenalty = gapExtensionPenalty
		self._mismatchPenalty = mismatchPenalty
		self._matchBonus = matchBonus
		self._matchMatrix = None
		self._gapMatrix1 = None
		self._gapMatrix2 = None
		self._bestAlignment = None

	def __repr__(self):
		"""
		Allow for a nice print representation of the optimal alignment.
		"""
		if self._bestAlignment is not None:
			return "s1: {0} \ns2: {1}".format(self._bestAlignment[0], self._bestAlignment[1])
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
		if self._matchMatrix is None:
			self._computeScoreMatrix()
		return min(self._matchMatrix[-1][-1], self._gapMatrix1[-1][-1], self._gapMatrix2[-1][-1])

	def _initializeMatrix(self):
		"""
		Initializes the score matrix for the alignment of s1 and s2

		s1 is placed on the rows of the matrix, while s2 is on the columns

		Helper function for computeScoreMatrix()
		"""
		# initialize an len(s1) + 1 x len(s2) + 1 matrix of 0s
		self._matchMatrix = zeroMatrix(len(self._s1) + 1, len(self._s2) + 1)
		self._gapMatrix1 = zeroMatrix(len(self._s1) + 1, len(self._s2) + 1)
		self._gapMatrix2 = zeroMatrix(len(self._s1) + 1, len(self._s2) + 1)

		# set the first row and first column of matchMatrix to -inf
		# set the first column of the gap matrices to -inf
		self._matchMatrix[0] = [float("inf") for j in range(len(self._s2) + 1)]
		# set the first row of the gap matrix for s1 to the gap penalty (including gap extension)
		self._gapMatrix1[0] = [self._gapOpeningPenalty + self._gapExtensionPenalty * (j - 1) for j in range(len(self._s2) + 1)]
		
		for i in range(len(self._s1) + 1):
			self._matchMatrix[i][0] = float("inf")
			self._gapMatrix1[i][0] = float("inf")
			self._gapMatrix2[i][0] = self._gapOpeningPenalty + self._gapExtensionPenalty * (i - 1)

		
		# gap matrix 2 is a transpose of gap matrix 1: gap penalty along column, undefined on row
		self._gapMatrix2[0] = [float("inf") for j in range(len(self._s2) + 1)]

		# print self._matchMatrix
		# print "\n"
		# print self._gapMatrix1
		# print "\n"
		# print self._gapMatrix2
		


	def _computeScoreMatrix(self):
		"""
		Fills in the score matrix for the alignment of s1 and s2

		Helper function for constructAlignment()
		"""
		self._initializeMatrix()
		for i in range(1, len(self._s1) + 1):
			for j in range(1, len(self._s2) + 1):
				
				if self._s1[i - 1] == self._s2[j - 1]:
					matchScore = self._matchBonus
				else:
					matchScore = self._mismatchPenalty

				self._matchMatrix[i][j] = matchScore + min(
					self._matchMatrix[i - 1][j - 1]
					, self._gapMatrix1[i - 1][j - 1]
					, self._gapMatrix2[i - 1][j - 1]
				)

				self._gapMatrix1[i][j] = min(
					self._gapOpeningPenalty + self._gapExtensionPenalty + self._matchMatrix[i][j - 1]
					, self._gapExtensionPenalty + self._gapMatrix1[i][j - 1]
					, self._gapOpeningPenalty + self._gapExtensionPenalty + self._gapMatrix2[i][j - 1]
				)

				self._gapMatrix2[i][j] = min(
					self._gapOpeningPenalty + self._gapExtensionPenalty + self._matchMatrix[i - 1][j]
					, self._gapOpeningPenalty + self._gapExtensionPenalty + self._gapMatrix1[i - 1][j]
					, self._gapExtensionPenalty + self._gapMatrix2[i - 1][j]
				)

		# print self._matchMatrix
		# print "\n"
		# print self._gapMatrix1
		# print "\n"
		# print self._gapMatrix2


	def _constructAlignment(self):
		"""
		Uses the completed score matrix to return the optimal alignment between s1 and s2

		Alignments are stored as a tuple of (alignment for s1, alignment for s2)
		"""
		if self._matchMatrix is None:
			self._computeScoreMatrix()

		i, j = len(self._s1), len(self._s2) # keep track of s1 and s2
		s1Align, s2Align = "", ""

		# determine which matrix to start in
		score = self.getAlignmentScore()
		if score == self._matchMatrix[-1][-1]:
			currentMatrix = self._matchMatrix
		elif score == self._gapMatrix1[-1][-1]:
			currentMatrix = self._gapMatrix1
		else:
			currentMatrix = self._gapMatrix2

		while i > 0 or j > 0:
			#print i, j, s1Align, s2Align, score

			if currentMatrix == self._matchMatrix:
				s1Align += self._s1[i - 1]
				s2Align += self._s2[j - 1]
				i -= 1
				j -= 1

			elif currentMatrix == self._gapMatrix1:
				s1Align += "-"
				s2Align += self._s2[j - 1]
				j -= 1

			elif currentMatrix == self._gapMatrix2:
				s2Align += "-"
				s1Align += self._s1[i - 1]
				i -= 1

			score = min(self._matchMatrix[i][j], self._gapMatrix1[i][j], self._gapMatrix2[i][j])
			if score == self._matchMatrix[i][j]:
				currentMatrix = self._matchMatrix
			elif score == self._gapMatrix1[i][j]:
				currentMatrix = self._gapMatrix1
			elif score == self._gapMatrix2[i][j]:
				currentMatrix = self._gapMatrix2
			else:
				raise ValueError('Current score does not match any matrix!')
			

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
	S1 = "CACATATTATTCACT"
	S2 = "CAGATTATTTCAT"
	aligner = StringAligner(S1, S2)
	aligner.align()
	print aligner
	#CACATATTATTCACT 
	#CAGATTATT-TCA-T
	
