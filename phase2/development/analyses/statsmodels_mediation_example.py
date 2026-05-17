import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.mediation import Mediation

# ---------------------------------------------------------
# 1. Generate Synthetic Data
# ---------------------------------------------------------
np.random.seed(42)
n_students = 200

# X (Independent Variable): Hours of study
study_hours = np.random.normal(loc=10, scale=2, size=n_students)

# M (Mediator): Practice problems completed
# We simulate M being heavily influenced by X
practice_problems = (
    5 + (2.5 * study_hours) + np.random.normal(loc=0, scale=3, size=n_students)
)

# Y (Dependent Variable): Exam score
# We simulate Y being heavily influenced by M, and slightly by X directly
exam_scores = (
    20
    + (0.8 * study_hours)
    + (1.2 * practice_problems)
    + np.random.normal(loc=0, scale=5, size=n_students)
)

# Combine into a DataFrame
df = pd.DataFrame(
    {
        "Study_Hours": study_hours,
        "Practice_Problems": practice_problems,
        "Exam_Score": exam_scores,
    }
)

# ---------------------------------------------------------
# 2. Build the Underlying Regression Models
# ---------------------------------------------------------
# Step A: The Mediator Model (X -> M)
# How does studying affect the number of practice problems done?
mediator_model = sm.OLS.from_formula("Practice_Problems ~ Study_Hours", data=df)

# Step B: The Outcome Model (X + M -> Y)
# How do BOTH studying and practice problems affect the exam score?
outcome_model = sm.OLS.from_formula(
    "Exam_Score ~ Study_Hours + Practice_Problems", data=df
)

# ---------------------------------------------------------
# 3. Run the Mediation Analysis
# ---------------------------------------------------------
# We pass the two models into the Mediation class.
# We must specify the independent variable ('Study_Hours')
# and the mediator ('Practice_Problems').
med = Mediation(
    outcome_model, mediator_model, exposure="Study_Hours", mediator="Practice_Problems"
)

# Fit the model using bootstrapping (default is 1000 simulations) to get standard errors
med_result = med.fit(n_rep=1000)

# Print the results
print(med_result.summary())
