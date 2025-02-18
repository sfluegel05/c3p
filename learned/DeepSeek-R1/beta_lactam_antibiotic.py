"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
</think>

Okay, I need to create a function called is_beta_lactam_antibiotic that checks if a SMILES string represents a beta-lactam antibiotic. The definition says it's an organonitrogen heterocyclic antibiotic containing a beta-lactam ring. 

First, the beta-lactam ring is a four-membered ring with one nitrogen and a carbonyl group. The SMARTS pattern for that should be something like a four-membered ring where one atom is nitrogen (like [#7]) and another is a carbonyl carbon (C=O). So the SMARTS could be "[#7]1-C-C(=O)-C-1" or similar. Wait, but the ring might have different connectivity. Maybe using a more flexible pattern like "[#7]1~*~C(=O)~*1" where ~ allows any bond type.

Next, the molecule must be an antibiotic. But how to check that programmatically? Maybe look for common antibiotic substructures or functional groups often present in antibiotics. Alternatively, since the problem states "organonitrogen heterocyclic", we can check for a heterocycle containing nitrogen. The beta-lactam ring itself is a heterocycle, so maybe that's sufficient if combined with other typical features.

So the steps would be:
1. Check for the presence of a beta-lactam ring.
2. Ensure it's an antibiotic. Since antibiotic activity is hard to determine from structure alone, maybe check for common substituents like a carboxylate group (common in penicillins and cephalosporins) or other groups often present in beta-lactams.

Wait, the examples given (like benzylpenicillin) have a carboxylate group. Maybe requiring a carboxylic acid or its derivative (like ester) could help. But some beta-lactams might have different groups. Alternatively, maybe just the presence of the beta-lactam ring is enough, given the problem's definition.

But the definition says "organonitrogen heterocyclic antibiotic". So the molecule must be an antibiotic (which is a biological activity), but how to determine that from structure? Perhaps we can't directly, so maybe the presence of the beta-lactam ring is the main criterion. However, some molecules might have a beta-lactam ring but not be antibiotics. But according to the problem statement, the class