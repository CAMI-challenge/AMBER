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
    ../amber.py -g input/gsa_mapping.binning input/goofy_hypatia_2 input/naughty_carson_2 -r input/unique_common.tsv -k "circular element" -p 1 -l "CONCOCT (CAMI), MaxBin 2.0.2 (CAMI)" -o output --stdout --silent
    """
    Then the exit code should be 0
    And the file "output/results.tsv" should exist
    And the stdout should contain:
    """
                    Tool  Accuracy (bp)  Accuracy (seq)  Adjusted Rand index (bp)  Adjusted Rand index (seq)  Average completeness (bp)  Average completeness (seq)  Average purity (bp)  Average purity (seq)  CAMI 1 F1 score (bp)  CAMI 1 F1 score (seq)  CAMI 1 average completeness (bp)  CAMI 1 average completeness (seq)  CAMI 1 std error of av. completeness (bp)  CAMI 1 std error of av. completeness (seq)  Completeness (bp)  Completeness (seq)  F1 score (bp)  F1 score (seq)  F1 score for sample (bp)  F1 score for sample (seq)  Misclassification rate (bp)  Misclassification rate (seq)  Percentage of binned bp  Percentage of binned sequences  Purity (bp)  Purity (seq)  Rand index (bp)  Rand index (seq)    Sample  Std error of av. completeness (bp)  Std error of av. completeness (seq)  Std error of av. purity (bp)  Std error of av. purity (seq) UniFrac (bp) UniFrac (seq)  avg_precision_bp_var  avg_recall_bp_var  avg_recall_bp_var_cami1 binning type rank
           Gold standard       1.000000        1.000000                  1.000000                   1.000000                   1.000000                    1.000000             1.000000              1.000000              1.000000               1.000000                          1.000000                           1.000000                                   0.000000                                    0.000000           1.000000            1.000000       1.000000        1.000000                   1.00000                   1.000000                     0.000000                      0.000000                 1.000000                        1.000000     1.000000      1.000000          1.00000          1.000000  CAMI_low                            0.000000                             0.000000                      0.000000                       0.000000         None          None              0.000000           0.000000                 0.000000       genome   NA
          CONCOCT (CAMI)       0.663234        0.321346                  0.643557                   0.750612                   0.933410                    0.742461             0.837135              0.757099              0.639181               0.544308                          0.516941                           0.424888                                   0.070242                                    0.056733           0.935822            0.381454       0.882655        0.749708                   0.79045                   0.503972                     0.315830                      0.257567                 0.966790                        0.432638     0.684170      0.742433          0.97174          0.945799  CAMI_low                            0.025427                             0.030531                      0.053293                       0.063328         None          None              0.073844           0.025862                 0.231892       genome   NA
     MaxBin 2.0.2 (CAMI)       0.809257        0.319959                  0.917149                   0.781649                   0.825796                    0.719064             0.947763              0.794919              0.866807               0.745751                          0.798593                           0.702312                                   0.058312                                    0.055000           0.837613            0.368405       0.882586        0.755092                   0.88359                   0.496961                     0.065092                      0.236671                 0.863569                        0.419098     0.934908      0.763329          0.99470          0.951167  CAMI_low                            0.050388                             0.050163                      0.016707                       0.044926         None          None              0.009211           0.101558                 0.136012       genome   NA

    """