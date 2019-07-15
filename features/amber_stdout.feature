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
                    Tool  >0.5compl<0.05cont  >0.5compl<0.1cont  >0.7compl<0.05cont  >0.7compl<0.1cont  >0.9compl<0.05cont  >0.9compl<0.1cont  Accuracy (bp)  Accuracy (seq)  Adjusted Rand index (bp)  Adjusted Rand index (seq)  Average completeness (bp)  Average completeness (seq)  Average purity (bp)  Average purity (seq)  Completeness (bp)  Completeness (seq)  Misclassification rate (bp)  Misclassification rate (seq)  Percentage of assigned bps  Percentage of assigned sequences  Purity (bp)  Purity (seq)  Rand index (bp)  Rand index (seq)    Sample  Std error of average completeness (bp)  Std error of average completeness (seq)  Std error of average purity (bp)  Std error of average purity (seq) UniFrac (bp) UniFrac (seq) binning type rank
          CONCOCT (CAMI)                  16                 18                  16                 17                  14                 15       0.663234        0.321346                  0.643540                   0.750520                   0.516941                    0.424888             0.837135              0.757099           0.935822            0.381454                     0.315830                      0.257567                    0.969401                          0.432638     0.684170      0.742433         0.971708          0.945658  CAMI_low                                0.070242                                 0.056733                          0.053293                           0.063328         None          None       genome   NA
     MaxBin 2.0.2 (CAMI)                  23                 28                  23                 28                  21                 24       0.809257        0.319959                  0.917148                   0.781552                   0.798593                    0.702312             0.947763              0.794919           0.837613            0.368405                     0.065092                      0.236671                    0.865600                          0.419098     0.934908      0.763329         0.994690          0.951013  CAMI_low                                0.058312                                 0.055000                          0.016707                           0.044926         None          None       genome   NA

    """