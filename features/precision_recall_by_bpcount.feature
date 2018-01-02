Feature: Run Scripts without error

  Background: Initial Setup
    Given I downloaded the scripts

  Scenario: Run precision_recall_per_bpcount.py

    Given I create the directory "input"
    And I copy the example data files
      | source                                     | dest        |
      |  test_query.binning     | input |
      |  gsa_mapping.binning      | input |
    When I run the command
    """
    ../src/precision_recall_by_bpcount.py -g input/gsa_mapping.binning input/test_query.binning
    """
    Then the exit code should be 0
    And the stdout should contain:
     """
     precision	recall
     0.700	0.839
"""
