# 1 HMSC Preparation

This folder contains scripts that assemble data and define the HMSC model objects used for Chapter 1.

## Scripts

- `S0_Data_Definitions.R`: loads and filters environmental data, occurrence data, traits, spatial coordinates, and study design objects; saves the preprocessed HMSC input data.
- `S1_Model_Definitions.R`: loads the prepared data, defines the HMSC model formulae and random levels, and writes unfitted model objects and testing metadata.

## Inputs

These scripts currently expect the original project data layout under `Data/`. During publication cleanup, these paths should be reviewed against the new `data-new/` structure and any private/raw inputs should be documented clearly.

## Outputs

The current scripts write preprocessed data and unfitted model objects into the original project layout, including `Data/` and `HmscOutputs/`. For the cleaned release, outputs should be routed to documented intermediate locations or described as externally archived if they are too large for GitHub.

Run these scripts from the repository root.
