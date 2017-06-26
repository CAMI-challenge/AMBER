Feature: Run Scripts without error

  Background: Initial Setup
    Given I downloaded the scripts

  Scenario: Run precision_recall_per_genome.py

    Given I download the file "https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz" to "gold_standard.fasta"
    And I create the directory "input"
    And I copy the example data files
      | source                                     | dest        |
      |  test_query.binning     | input |
      |  gsa_mapping.binning      | input |
    When I run the command
    """
    ../precision_recall_per_genome.py  -g input/gsa_mapping.binning -q input/test_query.binning -f gold_standard.fasta
    """
    Then the exit code should be 0
    And the stdout should contain:
     """
     1285_BH	1.0	1.0	4566533	4566533	4566533
     1030878	1.0	0.996440945639	3276530	3276530	3288233
     evo_1286_AP.008	0.223534141808	0.996316122005	15082564	3371468	3383934
    """