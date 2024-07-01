##  Computing-NCA-Parameters-for-Drugs

# Load the required library
library(PKNCA)

# Load the Theoph data from the datasets package
data(Theoph, package="datasets")

# Convert Theoph data to a data frame
theoph_df <- as.data.frame(Theoph)

# Display the first few rows of the Theoph data
print(head(theoph_df))

# Create the concentration object
conc_obj <- PKNCAconc(theoph_df, conc~Time|Subject)

# Extract unique dosing data (one row per dose)
d_dose <- unique(theoph_df[theoph_df$Time == 0, c("Dose", "Time", "Subject")])

# Display the dosing data
print(d_dose)

# Create the dose object
dose_obj <- PKNCAdose(d_dose, Dose~Time|Subject)

# Merge concentration and dose data automatically
data_obj_automatic <- PKNCAdata(conc_obj, dose_obj)
# Display PKNCA options
print(PKNCA.options("single.dose.aucs"))

# Display the intervals in the automatic data object
print(data_obj_automatic$intervals)

# Define manual intervals
intervals_manual <- data.frame(start=0, end=Inf, cmax=TRUE, tmax=TRUE, aucinf.obs=TRUE, auclast=TRUE)

# Merge concentration and dose data manually with specified intervals
data_obj_manual <- PKNCAdata(conc_obj, dose_obj, intervals=intervals_manual)

# Display the intervals in the manual data object
print(data_obj_manual$intervals)

# Compute parameters automatically
results_obj_automatic <- pk.nca(data_obj_automatic)

# Display the first few rows of the automatic results
print(head(as.data.frame(results_obj_automatic)))

# Summarize the automatic results
print(summary(results_obj_automatic))

# Compute parameters manually
results_obj_manual <- pk.nca(data_obj_manual)

# Display the first few rows of the manual results
print(head(as.data.frame(results_obj_manual)))

# Summarize the manual results
print(summary(results_obj_manual))

# Superposition for multiple dose PK
d_conc <- PKNCAconc(theoph_df, conc~Time|Subject)
conc_obj_multi <- PKNCAconc(
  superposition(d_conc, tau=168, dose.times=seq(0, 144, by=24), n.tau=1, check.blq=FALSE),
  conc~time|Subject
)

# Display the concentration object for multiple doses
print(conc_obj_multi)

# Create the dose object for multiple doses
dose_obj_multi <- PKNCAdose(expand.grid(Subject=unique(as.data.frame(conc_obj_multi)$Subject), time=seq(0, 144, by=24)), ~time|Subject)

# Display the dose object for multiple doses
print(dose_obj_multi)

# Create the data object for multiple doses
data_obj <- PKNCAdata(conc_obj_multi, dose_obj_multi)

# Display the intervals for multiple doses
print(data_obj$intervals[,c("Subject", "start", "end")])

# Define manual intervals for multiple doses
intervals_manual <- data.frame(start=c(0, 144), end=c(24, 168), cmax=TRUE, auclast=TRUE)

# Create the data object for multiple doses with manual intervals
data_obj <- PKNCAdata(conc_obj_multi, dose_obj_multi, intervals=intervals_manual)

# Display the intervals for manual multiple doses
print(data_obj$intervals)

# Compute parameters for multiple doses
results_obj <- pk.nca(data_obj)

# Display the results for multiple doses
print(results_obj)

# Summarize the results for multiple doses
print(summary(results_obj))


