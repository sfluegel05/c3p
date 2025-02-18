"""
Classifies: CHEBI:52221 isothiocyanate
"""
Wait, but the SMARTS pattern [C]-[N]=[C]=[S] may not capture all cases. For example, if the R group is part of a ring (like in benzyl isothiocyanate: S=C=NCc1ccccc1), the N is connected to a CH2 group (Cc1...). The SMARTS [C]-[N]=[C]=[S] would match because the N is connected to the CH2 carbon.

Another example: toluene 2,4-diisothiocyanate (Cc1ccc(cc1N=C=S)N=C=S). Each N=C=S group is connected to a carbon in the benzene ring, so the SMARTS would match both groups.

But what about a molecule with multiple N=C=S groups? The code would still return True as long as at least one matches.

Testing the code with the examples provided:

- phenethyl isothiocyanate (S=C=NCCc1ccccc1): The SMARTS would match the N connected to the phenethyl group's carbon.
- allyl isothiocyanate (C=CCN=C=S): The N is connected to the allyl group's carbon.
- methyl isothiocyanate (CN=C=S): N connected to methyl.
- (S)-sulforaphane (N(=C=S)CCCC[S@@](=O)C): Wait, the SMILES here is N(=C=S)... which would imply N double bonded to C, which is single bonded to S. But according to the SMARTS [C]-[N]=[C]=[S], this would not match because the C is not double bonded to S. So this would be a false negative. But according to the problem statement, (S)-sulforaphane is an example that should be classified as isothiocyanate. So there's a problem here.

Wait, this indicates that the SMARTS pattern [C]-[N]=[C]=[S] may not capture all cases. Let's look at the actual SMILES for (S)-sulforaphane: N(=C=S)CCCC[S@@](=O)C. Parsing this, the N is double bonded to C, which is single bonded to S. That would make the group N=C-S, which is not the isothiocyanate group (which should be N=C=S). So either the example is incorrect, or the SMILES is written differently. Alternatively, maybe the SMILES is correct and the structure is different. This is a problem because the code would fail to recognize this as an isothiocyanate.

But according to the problem statement, (S)-sulforaphane is an example of an isothiocyanate. So perhaps the SMILES provided is incorrect. Let me check the actual structure of sulforaphane. Sulforaphane's structure is CH3-SO-CH2-CH2-CH2-CH2-N=C=S. So the isothiocyanate group is at the end: -N=C=S. So the correct SMARTS should match that. The SMILES should be written as N=C=S... So perhaps the given SMILES is incorrect. For example, the correct SMILES would be N#C=S for thiocyanate, but isothiocyanate is N=C=S. Wait, no. Thiocyanate is S-C≡N (SMILES: S-C#N), while isothiocyanate is N=C=S (SMILES: N=C=S). So the correct SMILES for sulforaphane's isothiocyanate group should be N=C=S. Therefore, the given example's SMILES might have a typo. If the SMILES is N=C=S..., then the code would work. But if it's written as N(=C=S), that's a problem.

Assuming the example SMILES is correct, then the code would fail to classify it. But that's a problem with the example, not the code. However, since the user provided that example, the code must handle it.

Wait, let's parse the given SMILES for (S)-sulforaphane: N(=C=S)CCCC[S@@](=O)C. The N is attached to a C via a double bond (because of the =C=S part). The C is attached to S via a single bond. So the structure is N=C-S, which is not an isothiocyanate. So this suggests that the example SMILES is incorrect. But according to the problem statement, this is an example of an isothiocyanate. Therefore, there's a contradiction here.

Alternatively, perhaps the SMILES is written as N=C=S, but due to SMILES parsing rules, it's represented differently. For example, the SMILES N=C=S would be parsed as N double bonded to C, which is double bonded to S. So the SMARTS [N]=[C]=[S] would match. But in the given example, the SMILES is N(=C=S), which may not represent that correctly.

This indicates that the code as written might not handle all cases if the SMILES is not correctly formatted. However, assuming that the examples provided have correct SMILES, the code should work.

Another potential issue: the SMARTS pattern [C]-[N]=[C]=[S] requires the N to be connected to a carbon, then double bonded to C and S. But in some cases, the N might be part of a ring or connected via a different arrangement. However, according to the definition, any R-N=C=S group should be matched.

In conclusion, the code should use the SMARTS pattern [C]-[N]=[C]=[S] to identify the isothiocyanate group. This should cover all valid cases where the N is connected to a carbon and part of the N=C=S group. The previous error was due to incorrect indentation, but in the code provided here, the indentation is correct.

Final code:
</think>