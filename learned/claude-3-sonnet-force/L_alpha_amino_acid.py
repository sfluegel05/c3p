"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: CHEBI:32567 L-alpha-amino acid
Any alpha-amino acid having L-configuration at the alpha-carbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alpha-amino acid pattern (C-N-C(=O)-O)
    amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)C(=O)[O;H1,-]")
    matches = mol.GetSubstructMatches(amino_acid_pattern)

    if not matches:
        return False, "No alpha-amino acid pattern found"

    # Identify the alpha-carbon
    alpha_carbon = None
    for match in matches:
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == "C" and atom.GetHybridization() == Chem.HybridizationType.SP3:
            alpha_carbon = atom
            break

    if alpha_carbon is None:
        return False, "No alpha-carbon found"

    # Determine chirality at alpha-carbon
    try:
        conf = alpha_carbon.GetProp("_CIPCode")
        if conf == "R":
            return False, "D-configuration at alpha-carbon"
        elif conf == "S":
            return True, "L-configuration at alpha-carbon"
        else:
            # If CIP code is not available or unclear, analyze the 3D conformation (if available)
            if mol.GetConformer().Is3D():
                perm = list(alpha_carbon.GetProp("_CIPRank"))
                if perm == [1, 2, 3, 4]:
                    return True, "L-configuration at alpha-carbon (based on 3D conformation)"
                elif perm == [4, 3, 2, 1]:
                    return False, "D-configuration at alpha-carbon (based on 3D conformation)"
                else:
                    return False, "Unable to determine chirality at alpha-carbon (ambiguous 3D conformation)"
            else:
                return False, "Unable to determine chirality at alpha-carbon (no 3D conformation available)"
    except KeyError:
        # If CIP code and rank are not available, check for additional criteria
        if mol.GetConformer().Is3D():
            return False, "Unable to determine chirality at alpha-carbon (missing stereochemistry information)"
        else:
            # Check for additional criteria like molecular weight, atom counts, and functional groups
            mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
            if mol_wt < 50 or mol_wt > 400:
                return False, "Molecular weight outside typical range for alpha-amino acids"

            c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
            o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

            if c_count < 2 or n_count != 1 or o_count < 2:
                return False, "Atom counts inconsistent with alpha-amino acid structure"

            # Look for carboxyl and amino groups
            carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
            amino_pattern = Chem.MolFromSmarts("N")
            if not mol.HasSubstructMatch(carboxyl_pattern) or not mol.HasSubstructMatch(amino_pattern):
                return False, "Missing carboxyl or amino group"

            return True, "Meets basic structural criteria for alpha-amino acid (no stereochemistry information)"

# Improved error handling
try:
    result, reason = is_L_alpha_amino_acid("NC(CCO)C(O)=O")
    print(f"Result: {result}, Reason: {reason}")
except Exception as e:
    print(f"Error executing is_L_alpha_amino_acid('NC(CCO)C(O)=O'): {e}")