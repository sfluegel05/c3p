"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
Wait, but the code has errors. For example, `Chem.PathToSubmol` is not a valid function. Also, checking if a ring has oxygen by using a substructure match on the entire molecule may not be accurate. Let's correct that.

Revised code:

Instead of creating a submol, iterate over each atom in the ring to check for oxygen.

Also, the amino_pattern might need adjustment. Let's try: