"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python
"""
Classifies: Non-proteinogenic amino acid
Definition: Any amino-acid that is not naturally encoded in the genetic code of any organism.

This function examines whether a molecule (given by its SMILES string) qualifies as a free amino acid (containing at least one 
carboxyl group and at least one free amino group) and is not one of the 20 canonical proteinogenic amino acids.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Pre-compute InChIKeys for the 20 canonical proteinogenic amino acids.
# (The canonical amino acids are defined by their SMILES strings.)
_PROTEINOGENIC_SMILES = [
    "NCC(=O)O",                         # Glycine
    "CC(N)C(=O)O",                       # Alanine
    "CC(C)[C@H](N)C(=O)O",                # Valine
    "CCC[C@H](N)C(=O)O",                  # Leucine
    "CC[C@H](C)[C@H](N)C(=O)O",           # Isoleucine
    "C1CC[C@@H](NC1)C(=O)O",              # Proline
    "N[C@@H](Cc1ccccc1)C(=O)O",           # Phenylalanine
    "N[C@@H](Cc1ccc(O)cc1)C(=O)O",        # Tyrosine
    "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",   # Tryptophan
    "N[C@@H](CO)C(=O)O",                 # Serine
    "N[C@@H]([C@@H](O)C)C(=O)O",         # Threonine
    "N[C@@H](CS)C(=O)O",                 # Cysteine
    "N[C@@H](CCSC)C(=O)O",               # Methionine
    "N[C@@H](CC(=O)O)C(=O)O",            # Aspartic acid
    "N[C@@H](CCC(=O)O)C(=O)O",           # Glutamic acid
    "N[C@@H](CC(=O)N)C(=O)O",            # Asparagine
    "N[C@@H](CCC(=O)N)C(=O)O",           # Glutamine
    "N[C@@H](CCCCN)C(=O)O",              # Lysine
    "N[C@@H](CCCNC(=N)N)C(=O)O",         # Arginine
    "N[C@@H](Cc1cnc[nH]1)C(=O)O",         # Histidine
]
_PROTEINOGENIC_INCHIKEYS = set()
for smi in _PROTEINOGENIC_SMILES:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        # Note: We compute the InChIKey on the canonical form of the amino acid.
        ik = Chem.MolToInchiKey(mol)
        _PROTEINOGENIC_INCHIKEYS.add(ik)

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines whether the SMILES string represents a non-proteinogenic amino acid.
    A valid amino acid here must have:
      - At least one carboxyl group (either as a deprotonated carboxylate or a protonated carboxylic acid).
      - At least one free amino group that is not in an amide linkage.
      - An InChIKey that does not match one of the 20 canonical proteinogenic amino acids.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): (True, reason) if it qualifies as a non-proteinogenic amino acid;
                     (False, reason) otherwise.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens to reliably detect attached H atoms.
    mol = Chem.AddHs(mol)
    
    # Define two SMARTS patterns for a carboxyl group:
    # 1. Protonated carboxylic acid: C(=O)O (where O has at least one hydrogen)
    # 2. Deprotonated carboxylate: C(=O)[O-]
    carboxy_acid = Chem.MolFromSmarts("[CX3](=O)[O;H]")
    carboxy_deprot = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    # Get matches for either pattern.
    acid_matches = mol.GetSubstructMatches(carboxy_acid)
    deprot_matches = mol.GetSubstructMatches(carboxy_deprot)
    if not acid_matches and not deprot_matches:
        return False, "No carboxyl group found"
    
    # Search for at least one free amino group.
    # A free amino group is defined here as a nitrogen (atomic number 7) that has at least one hydrogen,
    # and which is not directly bonded to a carbonyl carbon.
    free_amino_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # Check that the nitrogen has at least one hydrogen.
        if atom.GetTotalNumHs() < 1:
            continue
        # Check if this nitrogen is directly bound to a carbon that is doubly-bonded to oxygen.
        is_amide_like = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # potential carbon
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                            is_amide_like = True
                            break
                if is_amide_like:
                    break
        if not is_amide_like:
            free_amino_found = True
            break
    if not free_amino_found:
        return False, "No free amino group found"
    
    # Compute the InChIKey for the molecule.
    try:
        inchi_key = Chem.MolToInchiKey(mol)
    except Exception as e:
        return False, f"Failed to generate InChIKey: {e}"
        
    # Check if the molecule matches one of the 20 canonical amino acids.
    if inchi_key in _PROTEINOGENIC_INCHIKEYS:
        return False, "Matches a canonical proteinogenic amino acid"
    
    return True, "Contains at least one carboxyl group and a free amino group and is non-proteinogenic"

# Example usage (to test the function, uncomment and run):
# test_smiles = [
#     "[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)O)=O)NC([C@@H](C3=CC=C(C=C3)OCC[C@H](C(O)=O)N)=O",  # nocardicin C
#     "NCCCC[C@H](N)CC(O)=O",  # (3S)-3,7-diaminoheptanoic acid
#     "[H][C@]1(O[C@@](C[C@H](O)[C@H]1NC(C)=O)(OC[C@H]1O[C@H](OC[C@H](N)C(O)=O)[C@H](NC(C)=O)[C@@H](O)[C@H]1O)C(O)=O)[C@H](O)[C@H](O)CO",  # O-[N-acetyl-alpha-neuraminyl-(2->6)-N-acetyl-alpha-D-galactosaminyl]-L-serine
# ]
# for smi in test_smiles:
#     result, reason = is_non_proteinogenic_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")