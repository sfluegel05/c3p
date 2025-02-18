"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
</think>

Okay, let's figure out why the previous code didn't work well. The main issues were false positives and false negatives. 

Looking at the false negatives first: Thymidine and others were missed because the deoxyribose detection failed. The original code looked for a 5-membered ring with one oxygen (like ribose), then checked if C2' had no hydroxyl. But maybe the SMARTS pattern for the sugar was too strict. Deoxyribose should have four carbons and one oxygen in the ring. Also, some structures might have substituents that confused the C2' check.

For the false positives, molecules like Scedapin C were incorrectly classified. The pyrimidine pattern "n1cnccc1" might match parts of larger rings that aren't actual pyrimidine bases. Also, the glycosidic bond check might have been too lenient, allowing connections to non-base nitrogens.

So, improvements needed:

1. **Better deoxyribose detection**: Use a SMARTS pattern that explicitly matches the deoxyribose structure (C1OCC(C)CO). Wait, no, deoxyribose is a furanose ring with four carbons and one oxygen. The correct SMARTS should be [C@H]1[C@H](O)C[C@@H](CO)O1 (but need to account for possible stereochemistry variations).

Wait, maybe using a more general deoxyribose pattern. The key is that in deoxyribose, the 2' carbon (C2) has only two substituents (H and the next carbon), whereas in ribose it has an OH. So the SMARTS for deoxyribose should be a furanose ring (O1C-C-C-C-C1) where the C2 has no oxygen attached.

2. **Pyrimidine base check**: The original pattern might miss substituted pyrimidines. A better approach is to look for a six-membered ring with two nitrogens, possibly with substitutions. A SMARTS like n1c([#6])ncc1 could be better, allowing substituents on the pyrimidine.

3. **Glycosidic bond verification**: Ensure the nitrogen in the pyrimidine is connected to the sugar via a glycosidic bond. The previous check might have allowed any connection, but it should specifically be between the base's nitrogen and the sugar's C1'.

4. **Phosphate check**: The original code checked for phosphorus, which is correct for nucleosides vs nucleotides.

So, revising the steps:

- Use a SMARTS pattern that accurately identifies the deoxyribose sugar, ensuring C2' has no hydroxyl.
- Use a more inclusive pyrimidine base pattern that allows various substituents.
- Strictly check that the glycosidic bond connects the base's nitrogen to the sugar's C1' carbon.
- Exclude molecules with phosphate groups.

Let me adjust the code accordingly. For the deoxyribose, perhaps a SMARTS like [C@H]1[C@H](O)C[C@@H](CO)O1 might not cover all cases due to stereochemistry variations. Alternatively, use a pattern that looks for the furanose ring with four carbons and one oxygen, and check that the C2' has no hydroxyl.

Also, the previous code's check for C2' hydroxyl was perhaps too strict. Maybe using a SMARTS that enforces the absence of an oxygen on C2'.

Putting this together, the improved code should:

1. Find a deoxyribose (furanose ring with O, four Cs, and C2' has no OH).
2. Find a pyrimidine base (six-membered ring with two nitrogens).
3. Ensure the base is connected via a glycosidic bond (N connected to C1' of the sugar).
4. No phosphate groups.

Testing this approach should reduce false positives and negatives.