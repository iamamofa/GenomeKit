from collections import Counter


class ConsensusBuilder:
    """
    Simple progressive Multiple Sequence Alignment (MSA)
    with consensus sequence generation.
    """

    def __init__(self, sequences_):

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

            bad = set(seq.upper()) - allowed

            if bad:
                raise ValueError(
                    f"Invalid characters found: {bad}"
                )

    def needleman_wunsch(self, seq1, seq2):
        """
        Global pairwise alignment using Needleman-Wunsch
        """

        match = 1
        mismatch = -1
        gap = -2

        n = len(seq1)
        m = len(seq2)

        # scoring matrix
        score = [
            [0] * (m + 1)
            for _ in range(n + 1)
        ]

        # traceback matrix
        trace = [
            [None] * (m + 1)
            for _ in range(n + 1)
        ]

        # initialize first column
        for i in range(1, n + 1):

            score[i][0] = i * gap
            trace[i][0] = "up"

        # initialize first row
        for j in range(1, m + 1):

            score[0][j] = j * gap
            trace[0][j] = "left"

        # fill matrices
        for i in range(1, n + 1):

            for j in range(1, m + 1):

                diag = score[i - 1][j - 1] + (
                    match
                    if seq1[i - 1] == seq2[j - 1]
                    else mismatch
                )

                up = score[i - 1][j] + gap

                left = score[i][j - 1] + gap

                best = max(diag, up, left)

                score[i][j] = best

                if best == diag:
                    trace[i][j] = "diag"

                elif best == up:
                    trace[i][j] = "up"

                else:
                    trace[i][j] = "left"

        # traceback
        aligned1 = []
        aligned2 = []

        i = n
        j = m

        while i > 0 or j > 0:

            if (
                i > 0
                and j > 0
                and trace[i][j] == "diag"
            ):

                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])

                i -= 1
                j -= 1

            elif (
                i > 0
                and (
                    j == 0
                    or trace[i][j] == "up"
                )
            ):

                aligned1.append(seq1[i - 1])
                aligned2.append("-")

                i -= 1

            else:

                aligned1.append("-")
                aligned2.append(seq2[j - 1])

                j -= 1

        return (
            "".join(reversed(aligned1)),
            "".join(reversed(aligned2))
        )

    def propagate_gaps(
        self,
        old_alignment,
        old_ref,
        new_ref
    ):
        """
        Insert new gaps into all previously aligned sequences
        based on updated reference alignment.
        """

        updated_alignment = []

        for seq in old_alignment:

            new_seq = []

            old_index = 0

            for char in new_ref:

                if char == "-":

                    new_seq.append("-")

                else:

                    new_seq.append(seq[old_index])
                    old_index += 1

            updated_alignment.append(
                "".join(new_seq)
            )

        return updated_alignment

    def build_msa(self):

        self._validate()

        # start with first sequence
        self.aligned = [self.sequences[0]]

        # progressively align remaining sequences
        for seq in self.sequences[1:]:

            current_ref = self.aligned[0]

            ref_aln, seq_aln = self.needleman_wunsch(
                current_ref,
                seq
            )

            # update previous aligned sequences
            updated_alignment = self.propagate_gaps(
                self.aligned,
                current_ref,
                ref_aln
            )

            # add new sequence
            updated_alignment.append(seq_aln)

            self.aligned = updated_alignment

        # reset profile
        self.profile = None

    def build_profile(self):

        if not self.aligned:
            self.build_msa()

        self.profile = []

        for column in zip(*self.aligned):

            counts = Counter(column)

            self.profile.append(counts)

    def consensus(self):

        if not self.aligned:
            self.build_msa()

        if self.profile is None:
            self.build_profile()

        consensus_seq = []

        for column in self.profile:

            counts = dict(column)

            # ignore gaps
            counts.pop("-", None)

            if not counts:

                consensus_seq.append("-")

            else:

                base = max(
                    counts,
                    key=counts.get
                )

                consensus_seq.append(base)

        return "".join(consensus_seq)

    def show_alignment(self):

        if not self.aligned:
            self.build_msa()

        for seq in self.aligned:
            print(seq)


# ---------------------------------------------------
# TESTS
# ---------------------------------------------------

print("\nTEST 1: Similar sequences")

seqs1 = [
    "ATAC",
    "ATAC",
    "ATGC"
]

cb1 = ConsensusBuilder(seqs1)

cb1.build_msa()

cb1.show_alignment()

print("Consensus:", cb1.consensus())


print("\nTEST 2: Small variations")

seqs2 = [
    "ATGGTA",
    "ATGCTA",
    "ATGCTA"
]

cb2 = ConsensusBuilder(seqs2)

cb2.build_msa()

cb2.show_alignment()

print("Consensus:", cb2.consensus())


print("\nTEST 3: Different lengths")

seqs3 = [
    "ATGGTA",
    "ATCT",
    "ATCCTA",
    "ATCA"
]

cb3 = ConsensusBuilder(seqs3)

cb3.build_msa()

cb3.show_alignment()

print("Consensus:", cb3.consensus())


print("\nTEST 4: Highly different lengths")

seqs4 = [
    "ATGCTAGCTA",
    "ATGC",
    "ATGCGCTA",
    "AT"
]

cb4 = ConsensusBuilder(seqs4)

cb4.build_msa()

cb4.show_alignment()

print("Consensus:", cb4.consensus())
