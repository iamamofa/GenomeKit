from collections import Counter

class ConsensusBuilder:
    def __init__(self,sequences_):
        self.sequences = sequences_
        self.aligned = []
        self.profile = None

    def _validate(self):
        if not self.sequences:
            raise ValueError("No sequences provided")

        if not all(isinstance(s, str) for s in self.sequences):
            raise TypeError("All sequences must be strings")

        allowed = set("ATGCN-")

        for seq in self.sequences:
            bad = set(seq) - allowed
            if bad:
                raise ValueError(f"Invalid characters found: {bad}")


    def needleman_wunsch(self, seq1, seq2):
        """
	    Performs a progressive Multiple Sequence Alignment for sequences
        """

        match, mismatch, gap = 1, -1, -2
        n,m = len(seq1), len(seq2)

        score =[[0]*(m+1) for _ in range(n+1)]
        trace = [[None]*(m+1) for _ in range(n+1)]

        for i in range(1, n+1):
            score[i][0] = i * gap
            trace[i][0] = "up"

        for j in range(1, m+1):
            score[0][j] = j * gap
            trace[0][j] = "left"

        for i in range(1, n+1):
            for j in range(1, m+1):
                diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
                up = score[i-1][j] + gap
                left = score[i][j-1] + gap

                best = max(diag, up, left)
                score[i][j] = best

                if best == diag:
                    trace[i][j] = "diag"
                elif best == up:
                    trace[i][j] = "up"
                else:
                    trace[i][j] = "left"

        a1, a2 = [], []
        i, j = n, m

        while i > 0 or j > 0:

            if i > 0 and j > 0 and trace[i][j] == "diag":
                a1.append(seq1[i-1])
                a2.append(seq2[j-1])
                i -= 1
                j -= 1

            elif i > 0 and (j == 0 or trace[i][j] == "up"):
                a1.append(seq1[i-1])
                a2.append("-")
                i -= 1

            else:
                a1.append("-")
                a2.append(seq2[j-1])
                j -= 1

        return "".join(reversed(a1)), "".join(reversed(a2))

    def build_msa(self):
    	        
        # start with first sequence
        #self.aligned = [self.sequences[0]]

        for seq in self.sequences[0:]:

            new_aligned = []

            #ref = self.aligned[0]
            ref = self.consensus() if self.aligned else self.sequences[0]

            ref_aln, seq_aln = self.needleman_wunsch(ref, seq)

            new_aligned.append(ref_aln)
            new_aligned.append(seq_aln)

            for old_seq in self.aligned[1:]:
                _, updated = self.needleman_wunsch(ref_aln, old_seq)
                new_aligned.append(updated)

            self.aligned = new_aligned

    def build_profile(self):

        self.profile = []

        for col in zip(*self.aligned, strict=False):
            counts = Counter(col)
            self.profile.append(counts)

    def consensus(self):
        self._validate()

        if not self.aligned:
            self.build_msa()

        if self.profile is None:
            self.build_profile()

        result = []

        for col in self.profile:
            col = dict(col)
            col.pop("-", None)

            if not col:
                result.append("-")
            else:
                result.append(max(col, key=col.get))

        return "".join(result)

seqs_identical = [
    "ATAC",
    "ATAC",
    "ATGC"
]

cb = ConsensusBuilder(seqs_identical)
print(cb.consensus())

seqs_variation = [
    "ATGGTA",
    "ATGCTA",
    "ATGCTA"
]
cb = ConsensusBuilder(seqs_variation)
print(cb.consensus())


seqs_diff_len = [
    "ATGGTA",
    "ATCT",
    "ATCCTA",
    "ATCA"
]

cb2 = ConsensusBuilder(seqs_diff_len)
print(cb2.consensus())