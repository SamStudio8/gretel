import unittest

EXP = [1,2,10]
def load_vcf(path, contig, start_pos, end_pos):
    from gretel import util
    return util.process_vcf(path, contig, start_pos, end_pos)

def load_bam(path, vcf_path, contig, start_pos, end_pos, threads):
    from gretel import util
    return util.load_from_bam(path, contig, start_pos, end_pos, load_vcf(vcf_path, contig, start_pos, end_pos), n_threads=threads)

class BasicRegressionTest(unittest.TestCase):

    def test_test(self):
        self.assertTrue(True)

    def test_vcf(self):
        VCF_h = load_vcf('tests/data/test.vcf.gz', 'hoot', 1, 19)
        self.assertEqual(VCF_h["N"], 3)
        self.assertEqual(VCF_h["snp_rev"], {0:1, 1:2, 2:10})
        self.assertEqual(VCF_h["snp_fwd"], {1:0, 2:1, 10:2})

        self.assertEqual(len(VCF_h["region"]), 19+1)
        for i in range(len(VCF_h["region"])):
            if i in EXP:
                self.assertEqual(VCF_h["region"][i], 1)
            else:
                self.assertEqual(VCF_h["region"][i], 0)
        return VCF_h

    def test_bam(self):

        for thread in [1,2]:
            hansel = load_bam('tests/data/test.bam', 'tests/data/test.vcf.gz', 'hoot', 1, 20, thread)

            self.assertEqual(hansel.n_slices, 5) # n reads
            self.assertEqual(hansel.n_crumbs, 9) # n snp connections
            self.assertTrue(hansel.L > 0) # Check L is actually set

            # Test a couple of elements in the Hansel matrix
            self.assertEqual(hansel.get_observation('_', 'A', 0, 1), 1)
            self.assertEqual(hansel.get_observation('A', 'A', 1, 2), 1)
            self.assertEqual(hansel.get_observation('A', 'A', 1, 3), 1)
            self.assertEqual(hansel.get_observation('A', 'A', 1, 4), 0)
            self.assertEqual(hansel.get_observation('C', 'C', 1, 2), 1)
            self.assertEqual(hansel.get_observation('C', 'C', 1, 3), 1)
            self.assertEqual(hansel.get_observation('C', 'C', 1, 4), 0)
            self.assertEqual(hansel.get_observation('T', 'T', 1, 2), 2)
            self.assertEqual(hansel.get_observation('G', 'G', 1, 2), 0)
            self.assertEqual(hansel.get_observation('G', 'G', 2, 3), 0)
            self.assertEqual(hansel.get_observation('G', 'G', 3, 4), 1)
            self.assertEqual(hansel.get_observation('G', '_', 4, 5), 1)


if __name__ == '__main__':
    unittest.main()
