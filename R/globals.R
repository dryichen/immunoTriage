# `triage_model` is a lazy-loaded data object in this package; declaring it
# here keeps R CMD check from reading the reference in triage_signature() as an
# undefined global.
utils::globalVariables("triage_model")
