"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python
"""
Classifies: Non-proteinogenic amino acid
Definition: Any amino-acid that is not naturally encoded in the genetic code of any organism.

The function is_non_proteinogenic_amino_acid takes a SMILES string as input and returns:
    (bool, str) -> (True, reason) if the molecule presents as a free amino acid (with at least one free carboxyl group
                     and at least one free amino group) and is not one of the 20 canonical proteinogenic amino acids;
                     (False, reason) otherwise.

Revised Heuristics:
  1. Parse the SMILES and add explicit hydrogens.
  2. Look for a carboxyl group using a SMARTS pattern that allows either –OH or a deprotonated –O⁻.
  3. Look for a free amino group defined as a nitrogen atom with at least one hydrogen that is not directly bound
     to a carbonyl carbon (i.e. avoiding an amide-bond situation).
  4. Compute the InChIKey and reject if the molecule matches one of the 20 canonical amino acids.
  
Note: Some non-proteinogenic amino acids contain additional amide bonds (intramolecular or in cyclic structures)
that are inherent to their backbone. Therefore, we do not outright reject the molecule for containing any amide linkages;
we simply require that the free acid and free amine functionalities can be found.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Pre-compute InChIKeys for the 20 canonical proteinogenic amino acids.
_PROTEINOGENIC_SMILES = [
    "NCC(=O)O",                         # Glycine
    "CC(N)C(=O)O",                       # Alanine
    "CC(C)[C@H](N)C(=O)O",                # Valine
    "CCC[C@H](N)C(=O)O",                  # Leucine
    "CC[C@H](C)[C@H](N)C(=O)O",           # Isoleucine
    "C1CC[C@@H](NC1)C(=O)O",              # Proline
    "N[C@@H](Cc1ccccc1)C(=O)O",           # Phenylalanine
    "N[C@@H](Cc1ccc(O)cc1)C(=O)O",         # Tyrosine
    "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",    # Tryptophan
    "N[C@@H](CO)C(=O)O",                  # Serine
    "N[C@@H]([C@@H](O)C)C(=O)O",          # Threonine
    "N[C@@H](CS)C(=O)O",                  # Cysteine
    "N[C@@H](CCSC)C(=O)O",                # Methionine
    "N[C@@H](CC(=O)O)C(=O)O",             # Aspartic acid
    "N[C@@H](CCC(=O)O)C(=O)O",            # Glutamic acid
    "N[C@@H](CC(=O)N)C(=O)O",             # Asparagine
    "N[C@@H](CCC(=O)N)C(=O)O",            # Glutamine
    "N[C@@H](CCCCN)C(=O)O",               # Lysine
    "N[C@@H](CCCNC(=N)N)C(=O)O",          # Arginine
    "N[C@@H](Cc1cnc[nH]1)C(=O)O",          # Histidine
]
_PROTEINOGENIC_INCHIKEYS = set()
for smi in _PROTEINOGENIC_SMILES:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        ik = Chem.MolToInchiKey(mol)
        _PROTEINOGENIC_INCHIKEYS.add(ik)

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines whether the SMILES string represents a free amino acid that is
    non-proteinogenic (i.e. not one of the 20 canonical amino acids), using revised heuristics.
    
    The function checks:
      - That the molecule parses and explicit hydrogens are added.
      - That the molecule contains at least one carboxyl group. Here we use a SMARTS pattern that
        allows matching either a protonated carboxylic acid or a deprotonated carboxylate.
      - That the molecule contains at least one free amino group. A free amino group is defined as a
        nitrogen (atomic number 7) with at least one hydrogen, where the nitrogen is not directly bonded
        to a carbonyl carbon (to avoid classifying an amide as a free amine).
      - That the molecule’s InChIKey is not found among the 20 canonical proteinogenic amino acids.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): (True, reason) if classified as a non‐proteinogenic amino acid;
                     (False, reason) otherwise.
    """
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens.
    mol = Chem.AddHs(mol)
    
    # Look for at least one carboxyl group.
    # The SMARTS below matches a carboxyl functional group including both protonated ([OX1]) and deprotonated ([O-]) forms.
    carboxyl_smarts = "[CX3](=O)[OX1,$([O-])]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"
    
    # Look for at least one free amino group.
    # We will loop over nitrogen atoms and require that:
    #   - They have at least one attached hydrogen.
    #   - They are not directly bonded to a carbonyl carbon.
    free_amino_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # Check number of attached hydrogens; use explicit+implicit count.
        if atom.GetTotalNumHs() < 1:
            continue
        # Check if this nitrogen is bonded to any carbon that is doubly-bonded to oxygen.
        is_amide_like = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon candidate
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
    
    # Compute the InChIKey of the hydrogen-added molecule.
    try:
        inchi_key = Chem.MolToInchiKey(mol)
    except Exception as e:
        return False, f"Failed to generate InChIKey: {e}"
    
    if inchi_key in _PROTEINOGENIC_INCHIKEYS:
        return False, "Matches a canonical proteinogenic amino acid"
    
    return True, "Contains at least one carboxyl group and a free amino group and is non-proteinogenic"

# Example usage (uncomment to test):
# test_smiles = [
#     "[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)O)=O)NC([C@@H](C3=CC=C(C=C3)OCC[C@H](C(O)=O)N)=O",  # nocardicin C 
#     "NCCCC[C@H](N)CC(O)=O",  # (3S)-3,7-diaminoheptanoic acid 
#     "[H][C@]1(O[C@@](C[C@H](O)[C@H]1NC(C)=O)(OC[C@H]1O[C@H](OC[C@H](N)C(O)=O)[C@H](NC(C)=O)[C@@H](O)[C@H]1O)C(O)=O)[C@H](O)[C@H](O)CO",  # O-[N-acetyl-alpha-neuraminyl-(2->6)-N-acetyl-alpha-D-galactosaminyl]-L-serine 
#     "N[C@@H](CC1=CC(=O)C(O)=CC1=O)C(O)=O",  # L-topaquinone 
#     "C[C@H](N)C[C@@H](N)C(O)=O",  # (2R,4S)-2,4-diaminopentanoic acid 
# ]
# for smi in test_smiles:
#     result, reason = is_non_proteinogenic_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")