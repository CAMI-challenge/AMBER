Feature: Run Scripts without error

  Background: Initial Setup
    Given I downloaded the scripts

  Scenario: Run html_plots.py
    Given I create the directory "input"
    And I copy the example data files
      | source                                     | dest        |
      |  precision_recall.tsv     | input |
      |  summary.tsv      | input |
    When I run the command
    """
    ../src/html_plots.py  -n "goofy" -o "out.html" -p "input/precision_recall.tsv" -s "input/summary.tsv"
    """
    Then the exit code should be 0
    And the file "out.html" should exist
