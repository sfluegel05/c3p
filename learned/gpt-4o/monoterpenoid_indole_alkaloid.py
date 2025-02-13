"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
Based on the failures observed in the previous attempt, it seems that the criteria for detecting monoterpenoid indole alkaloids (MIAs) were not comprehensive enough. The key issues were:

1. **Indole Moiety Detection**: The program failed to detect the indole moiety in several compounds even though they belong to the MIA class. It's crucial to ensure the SMARTS pattern for identifying the indole is correctly defined and covers the structural variations of the indole present in these molecules.

