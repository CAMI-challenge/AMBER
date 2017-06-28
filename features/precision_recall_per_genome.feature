Feature: Run Scripts without error

  Background: Initial Setup
    Given I downloaded the scripts

  Scenario: Run precision_recall_per_genome.py

    Given I download the file "https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz" to "gold_standard.fasta.gz"
    And I create the directory "input"
    And I copy the example data files
      | source                                     | dest        |
      |  test_query.binning     | input |
      |  gsa_mapping.binning      | input |
    When I run the command
    """
    ../precision_recall_per_genome.py -g input/gsa_mapping.binning -f gold_standard.fasta.gz input/test_query.binning
    """
    Then the exit code should be 0
    And the stdout should contain:
    """
    1285_BH	1.0	1.0	4566533	4566533	4566533
"""