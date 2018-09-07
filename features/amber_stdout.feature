Feature: Run Scripts without error

  Background: Initial Setup
    Given I downloaded the scripts

  Scenario: Run amber.py

    Given I create the directory "input"
    And I copy the example data files
      | source                                     | dest        |
      |  gsa_mapping.binning     | input |
      |  goofy_hypatia_2      | input |
      |  naughty_carson_2      | input |
      |  unique_common.tsv      | input |
    When I run the command
    """
    ../amber.py -g input/gsa_mapping.binning input/goofy_hypatia_2 input/naughty_carson_2 -r input/unique_common.tsv -k "circular element" -p 1 -l "CONCOCT (CAMI), MaxBin 2.0.2 (CAMI)" -o output --stdout
    """
    Then the exit code should be 0
    And the file "output/results.tsv" should exist
    And the stdout should contain:
    """
    tool binning type rank  average purity  standard deviation of average purity  standard error of average purity  average completeness  standard deviation of average completeness  standard error of average completeness  average purity per bp  average completeness per bp  accuracy  percentage of assigned bps  Rand index by bp counts  Rand index by sequence counts  adjusted Rand index by bp counts  adjusted Rand index by sequence counts  misclassification rate  >0.5compl<0.1cont  >0.7compl<0.1cont  >0.9compl<0.1cont  >0.5compl<0.05cont  >0.7compl<0.05cont  >0.9compl<0.05cont
         CONCOCT (CAMI)       genome   NA        0.837135                              0.271742                          0.053293              0.516941                                    0.481552                                0.070242               0.684170                     0.935822  0.663234                    0.969401                 0.971708                       0.945658                          0.643540                                0.750520                0.315830                 18                 17                 15                  16                  16                  14
    MaxBin 2.0.2 (CAMI)       genome   NA        0.947763                              0.095974                          0.016707              0.798593                                    0.368798                                0.058312               0.934908                     0.837613  0.809257                    0.865600                 0.994690                       0.951013                          0.917148                                0.781552                0.065092                 28                 28                 24                  23                  23                  21

    """