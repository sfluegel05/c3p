"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python
"""
Classifies: Non-proteinogenic amino acid
Definition: Any amino-acid that is not naturally encoded in the genetic code of any organism.
This function checks that the molecule contains at least one carboxyl group (protonated or deprotonated)
and at least one free amino group (i.e. an amine that is not in an amide linkage) and further that the molecule
does not match one of the 20 canonical proteinogenic amino acids.
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
        # Compute InChIKey on the canonical form (without extra hydrogens)
        canon = Chem.RemoveHs(mol)
        ik = Chem.MolToInchiKey(canon)
        _PROTEINOGENIC_INCHIKEYS.add(ik)

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines whether the SMILES string represents a non-proteinogenic amino acid.
    A valid amino acid must have:
      - At least one carboxyl group (either as a protonated carboxylic acid or a deprotonated carboxylate)
      - At least one free amino group (an amine group not involved in an amide linkage)
      - An InChIKey that does not match one of the 20 canonical proteinogenic amino acids.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): (True, reason) if the molecule qualifies as a non-proteinogenic amino acid;
                     (False, reason) otherwise.
    """
    # Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens for reliable substructure detection.
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns for a carboxyl group:
    # Protonated: e.g. C(=O)O with at least one attached hydrogen
    carboxy_acid = Chem.MolFromSmarts("[CX3](=O)[O;H]")
    # Deprotonated: e.g. C(=O)[O-]
    carboxy_deprot = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    acid_matches = mol.GetSubstructMatches(carboxy_acid)
    deprot_matches = mol.GetSubstructMatches(carboxy_deprot)
    if not acid_matches and not deprot_matches:
        return False, "No carboxyl group found"
    
    # Define a SMARTS pattern for a 'free' amino group.
    # This pattern matches a trivalent nitrogen with one or two hydrogens that is not directly bonded to a carbonyl carbon.
    free_amino_pattern = Chem.MolFromSmarts("[NX3;H1,H2;!$(N[C;$(C=O)])]")
    if not mol.HasSubstructMatch(free_amino_pattern):
        return False, "No free amino group found"
    
    # Remove explicit hydrogens for canonical InChIKey generation.
    canon_mol = Chem.RemoveHs(mol)
    
    try:
        inchi_key = Chem.MolToInchiKey(canon_mol)
    except Exception as e:
        return False, f"Failed to generate InChIKey: {e}"
    
    if inchi_key in _PROTEINOGENIC_INCHIKEYS:
        return False, "Matches a canonical proteinogenic amino acid"
    
    return True, "Contains at least one carboxyl group and a free amino group and is non-proteinogenic"

# Example usage:
# Uncomment the following lines to test with a few sample compounds.
# test_smiles = [
#     "[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)O)=O)NC([C@@H](C3=CC=C(C=C3)OCC[C@H](C(O)=O)N)=O",  # nocardicin C
#     "NCCCC[C@H](N)CC(O)=O",  # (3S)-3,7-diaminoheptanoic acid
#     "N[C@@H](CC1=CC(=O)C(O)=CC1=O)C(O)=O",  # L-topaquinone
# ]
# for smi in test_smiles:
#     result, reason = is_non_proteinogenic_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")