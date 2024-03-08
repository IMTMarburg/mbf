import mbf.genomes
import tempfile
import pytest
import mbf.functional.databases as databases
import pypipegraph2 as ppg


class TestsGroupParsing:
    def do_check(self, content, supposed, genome):
        tf = tempfile.NamedTemporaryFile(mode='w')
        tf.write(content)
        tf.flush()
        parser = databases.GroupsFromFile("tempfile", tf.name)
        actual = parser.get_sets(genome)
        assert actual == supposed

    def test_hugo(self,use_prebuild_genome):
        human_genome = mbf.genomes.EnsemblGenome("Homo_sapiens", 108)
        file = """
        #I am ignored
        !HUGO testgroup
        PPARD, PPARA
        PPARG
        !HUGO testgroup 2
        Angptl4 PPARD
        """
        supposed = {
            "testgroup": set(
                [
                    "ENSG00000112033",
                    "ENSG00000186951",
                    "ENSG00000132170",
                ]
            ),
            "testgroup 2": set(["ENSG00000167772", "ENSG00000112033"]),
        }
        self.do_check(file, supposed, human_genome)

    def test_hugo_raises(self, use_prebuild_genome):
        human_genome = mbf.genomes.EnsemblGenome("Homo_sapiens", 108)
        file = """
        #I am ignored
        !HUGO testgroup
        PPARD, PPARA
        PPARG
        !HUGO testgroup 2
        Angptl4 PPARD doesnotexist
        """
        supposed = {
            "testgroup": set(
                [
                    "ENSG00000112033",
                    "ENSG00000186951",
                    "ENSG00000132170",
                ]
            ),
            "testgroup 2": set(["ENSG00000167772", "ENSG00000112033"]),
        }

        with pytest.raises(ValueError):
            self.do_check(file, supposed, human_genome)

    def test_stable_id_human(self, use_prebuild_genome):
        human_genome = mbf.genomes.EnsemblGenome("Homo_sapiens", 108)
        file = """
        !stable_id testgroup
        ENSG00000112033, ENSG00000186951
        ENSG00000132170
        """
        supposed = {
            "testgroup": set(
                [
                    "ENSG00000112033",
                    "ENSG00000186951",
                    "ENSG00000132170",
                ]
            ),
        }
        self.do_check(file, supposed, human_genome)

    @pytest.mark.xfail
    def test_stable_id_human_translated_to_mouse(self, use_prebuild_genome):
        mouse_genome = mbf.genomes.EnsemblGenome("Mus_musculus", 108)
        file = """
        !stable_id testgroup
        ENSG00000112033, ENSG00000186951
        ENSG00000132170
        """
        supposed = {
            "testgroup": set(
                [
                    "ENSMUSG00000000440",
                    "ENSMUSG00000002250",
                    "ENSMUSG00000022383",
                ]
            ),
        }
        self.do_check(file, supposed, mouse_genome)

    def test_stable_id_mouse(self, use_prebuild_genome):
        mouse_genome = mbf.genomes.EnsemblGenome("Mus_musculus", 108)
        file = """
        !stable_id testgroup
        ENSMUSG00000000440,ENSMUSG00000002250,ENSMUSG00000022383
        """
        supposed = {
            "testgroup": set(
                [
                    "ENSMUSG00000000440",
                    "ENSMUSG00000002250",
                    "ENSMUSG00000022383",
                ]
            ),
        }
        self.do_check(file, supposed, mouse_genome)

    def test_mouse_name(self, use_prebuild_genome):
        mouse_genome = mbf.genomes.EnsemblGenome("Mus_musculus", 108)
        file = """
        !mouse_name testgroup
        PPARD, PPARA, PPARG
        """
        supposed = {
            "testgroup": set(
                [
                    "ENSMUSG00000000440",
                    "ENSMUSG00000002250",
                    "ENSMUSG00000022383",
                ]
            ),
        }
        self.do_check(file, supposed, mouse_genome)
