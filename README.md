# stringalign
Package to find the optimal alignment of two strings, s1 and s2

## Example

For the input strings ```s1 = AGGCT``` and ```s2 = AGCA```, the optimal alignment is
```("AGGCT", "AG-CA")```.

## Use

To align s1 and s2, you can directly use the ```StringAligner``` class:
```
> aligner = StringAligner(s1, s2)
> aligner.align() # returns tuple (best alignment for s1, best alignment for s2)
```
Or you can use the provided helper function align(s1, s2), which is a thin wrapper on StringAligner

If ```.align()``` has been called, you can call ```print``` on your aligner object to see a formatted version of the optimal alignment.

The alignment score for the optimal alignment is accessible via the ```.getAlignmentScore()``` method.


StringAligner takes two optional named parameters, mismatchPenalty and gapPenalty. When creating the object you can override the default values, but since the algorithm finds the best alignment through minimization, they must be strictly greater than zero. mismatchPenalty is the penalty for two unequal characters, while gapPenalty is the penalty for placing a gap in one string.

## Implementation

```stringalign``` uses the [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) using the [Affine Gap extension](https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf). Runtime is quadratic; O( max(len(s1), len(s2)) ** 2 ).
