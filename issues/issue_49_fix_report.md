# Issue #49 Fix Report: Enhanced Multicollinearity Detection in add_integration()

## Issue Summary

**Issue**: [#49](https://github.com/dmphillippo/multinma/issues/49) - Improve add_integration function to detect and report multicollinearity issues

**Problem**: The `add_integration()` function in multinma performs quasi-Monte Carlo numerical integration for multilevel network meta-regression (ML-NMR). While it had basic validation for correlation matrices (checking positive definiteness), it lacked comprehensive multicollinearity detection that could affect integration quality and model performance. Users needed better warnings when:

- Covariates are highly correlated (potential redundancy)
- Correlation matrices are near-singular (numerical instability)
- Multiple multicollinearity issues exist simultaneously
- Integration quality might be compromised by correlation structure

## Solution Implementation

### Enhanced Multicollinearity Detection

#### New Function: `detect_multicollinearity()`

Added a comprehensive internal function that detects multiple types of multicollinearity issues:

**Key Features:**
1. **High Pairwise Correlations**: Detects correlations above threshold (default |r| > 0.95)
2. **Condition Number Analysis**: Calculates ratio of largest to smallest eigenvalue (default threshold > 1000)
3. **Near-Singularity Detection**: Identifies correlation matrices with very small eigenvalues (< 1e-10)
4. **Determinant Analysis**: Additional check for matrix singularity using determinant
5. **Informative Reporting**: Specific variable names and correlation values in warnings
6. **Remedial Suggestions**: Actionable advice for addressing detected issues

#### Function Signature:
```r
detect_multicollinearity(cor_matrix, var_names, 
                        threshold_high = 0.95, 
                        threshold_condition = 1000)
```

### Integration Points for Detection

The enhanced multicollinearity detection is now integrated at **three critical points** in the `add_integration()` workflow:

1. **User-Provided Correlation Matrix**: Direct validation when users specify `cor` argument
2. **IPD-Derived Correlation Matrix**: Validation of weighted correlation matrix computed from IPD studies  
3. **Adjusted Correlation Matrix**: Validation after copula correlation adjustments

### Warning Types and Messages

#### 1. High Pairwise Correlations
```
Warning: High correlations detected (|r| > 0.95):
x1 & x2 (r = 0.97), x3 & x4 (r = 0.98)
Consider removing redundant covariates or using dimension reduction techniques.
```

#### 2. High Condition Number
```
Warning: Correlation matrix has high condition number (15432.1 > 1000).
This indicates potential multicollinearity issues that may affect integration quality.
Consider: (1) removing highly correlated covariates, (2) using principal components, 
or (3) checking data quality.
```

#### 3. Near-Singular Matrix
```
Warning: Correlation matrix is near-singular (smallest eigenvalue = 1.23e-12).
This may cause numerical instability in integration.
Check for perfectly correlated or redundant covariates.
```

#### 4. Low Determinant
```
Note: Correlation matrix determinant is very small (2.45e-11), indicating potential redundancy.
```

#### 5. Well-Conditioned Matrix
```
Info: Correlation matrix appears well-conditioned for numerical integration.
```

### Files Modified

#### Enhanced: `R/integration.R` (88 new lines)

**New Function:**
- `detect_multicollinearity()`: Comprehensive multicollinearity detection (76 lines)
- `plot_integration_error()`: Placeholder for integration diagnostic plots (12 lines)

**Enhanced Validation:**
- Added multicollinearity checks to `add_integration.data.frame()` method
- Added multicollinearity checks to `add_integration.nma_data()` method
- Integrated detection at all correlation matrix creation/validation points

#### New Tests: `tests/testthat/test-issue49-multicollinearity-detection.R` (263 lines)

**Test Coverage:**
1. High pairwise correlations detection
2. Near-singular matrix detection  
3. Perfect correlation detection
4. Well-conditioned matrix positive feedback
5. IPD-derived correlation validation
6. Multiple simultaneous issues
7. Direct function testing
8. Custom threshold functionality

### Technical Implementation Details

#### Multicollinearity Detection Algorithm

1. **Pairwise Correlation Analysis**:
   ```r
   upper_tri <- upper.tri(cor_matrix)
   high_corr_indices <- which(abs(cor_matrix[upper_tri]) > threshold_high)
   ```

2. **Condition Number Calculation**:
   ```r
   eigenvals <- eigen(cor_matrix, symmetric = TRUE, only.values = TRUE)$values
   condition_number <- max(eigenvals) / min(eigenvals)
   ```

3. **Singularity Detection**:
   ```r
   min_eigenval <- min(eigenvals)
   det_val <- det(cor_matrix)
   ```

#### Integration with Existing Workflow

The multicollinearity detection seamlessly integrates with the existing `add_integration()` workflow without breaking backward compatibility:

- Maintains all existing functionality
- Adds warnings/messages without changing return values
- Provides actionable feedback for users
- Handles edge cases (single variables, etc.)

### Benefits for Users

#### 1. **Early Problem Detection**
Users receive immediate feedback about correlation structure issues before costly model fitting.

#### 2. **Specific Variable Identification** 
Warnings identify exactly which variable pairs are problematic:
```
High correlations detected: age & age_squared (r = 0.987)
```

#### 3. **Actionable Guidance**
Warnings include specific suggestions for remediation:
- Remove redundant covariates
- Use principal component analysis
- Check data quality
- Consider dimension reduction

#### 4. **Integration Quality Assurance**
Helps ensure numerical integration accuracy by identifying issues that could affect QMC integration.

#### 5. **Model Performance Optimization**
Reduces likelihood of convergence issues and improves model reliability.

### Testing and Validation

#### Comprehensive Test Suite
- **8 test cases** covering all detection scenarios
- **Edge case testing**: Single variables, perfect correlations
- **Integration testing**: Real workflow with IPD and AgD data
- **Custom threshold testing**: Verifies flexible thresholds work
- **Multiple issue detection**: Simultaneous problems handled correctly

#### Example Test Cases
```r
# High correlation detection
expect_warning(add_integration(...), "High correlations detected.*x1 & x2.*r = 0.97")

# Condition number detection  
expect_warning(add_integration(...), "high condition number")

# Well-conditioned matrix feedback
expect_message(add_integration(...), "well-conditioned")
```

### Backward Compatibility

- **No breaking changes**: All existing code continues to work
- **Optional warnings**: Enhanced feedback without changing function signatures
- **Existing validation preserved**: All original correlation matrix checks remain
- **Graceful degradation**: Functions properly even with problematic matrices

### Usage Examples

#### Basic Usage (Automatic Detection)
```r
network <- add_integration(network,
                          age = distr(qgamma, mean = age_mean, sd = age_sd),
                          bmi = distr(qnorm, mean = bmi_mean, sd = bmi_sd),
                          smoking = distr(qbern, prob = smoking_prop))
# Automatically detects multicollinearity in IPD-derived correlations
```

#### With User-Provided Correlations
```r
cor_matrix <- matrix(c(1, 0.97, 0.97, 1), nrow = 2)  # High correlation
network <- add_integration(network,
                          age = distr(qgamma, mean = age_mean, sd = age_sd), 
                          age_squared = distr(qgamma, mean = age_sq_mean, sd = age_sq_sd),
                          cor = cor_matrix)
# Warning: High correlations detected (|r| > 0.95): age & age_squared (r = 0.97)
```

### Future Enhancements

The implementation provides a foundation for additional integration quality features:

1. **Integration Error Plotting**: Framework for `plot_integration_error()` function
2. **Variance Inflation Factors**: Could add VIF calculations for regression contexts
3. **Principal Component Suggestions**: Could automatically suggest PC transformations
4. **Adaptive Integration**: Could adjust integration points based on correlation structure

## Conclusion

This enhancement significantly improves the robustness and user-friendliness of the `add_integration()` function by providing comprehensive multicollinearity detection and informative feedback. Users can now:

- Identify problematic correlation structures before model fitting
- Receive specific guidance on which variables are causing issues  
- Take proactive steps to improve integration quality and model performance
- Benefit from better numerical stability in ML-NMR analyses

The implementation maintains full backward compatibility while providing valuable new functionality that enhances the overall quality of multinma analyses. 