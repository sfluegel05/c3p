from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonyl_groups(smiles: str):
    """
    Determines if a molecule contains a sulfonyl group (SO2 group).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains sulfonyl group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'S']
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    for sulfur in sulfur_atoms:
        # Get neighbors of sulfur
        neighbors = sulfur.GetNeighbors()
        
        # Count double-bonded oxygens
        double_bonded_oxygens = 0
        for neighbor in neighbors:
            if neighbor.GetSymbol() == 'O':
                # Get bond between S and O
                bond = mol.GetBondBetweenAtoms(sulfur.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_bonded_oxygens += 1

        if double_bonded_oxygens == 2:
            # Get substituents
            substituents = []
            for neighbor in neighbors:
                if neighbor.GetSymbol() != 'O':
                    if neighbor.GetSymbol() == 'C':
                        # Check if it's part of an aromatic ring
                        if neighbor.GetIsAromatic():
                            substituents.append("aryl")
                        else:
                            substituents.append("alkyl")
                    else:
                        substituents.append(neighbor.GetSymbol())
            
            # Classify specific types of sulfonyl groups
            if "aryl" in substituents:
                # Check for specific aromatic substituents
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1:
                        return True, "Nosyl group (nitrobenzenesulfonyl)"
                    elif atom.GetSymbol() == 'Br':
                        return True, "Brosyl group (p-bromobenzenesulfonyl)"
                    elif atom.GetSymbol() == 'N' and len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']) == 3:
                        return True, "Dansyl group (dimethylaminonaphthalenesulfonyl)"
                
                if len(substituents) == 2 and substituents.count("aryl") == 1:
                    if any(atom.GetSymbol() == 'C' and atom.GetIsAromatic() for atom in mol.GetAtoms()):
                        return True, "Tosyl group (p-toluenesulfonyl)"
                    return True, "Phenylsulfonyl group"
            
            return True, "Generic sulfonyl group"

    return False, "No sulfonyl groups (SO2) found"
# Pr=1.0
# Recall=1.0