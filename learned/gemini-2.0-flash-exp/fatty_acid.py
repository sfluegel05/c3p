"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import GraphDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    Fatty acids are characterized by a long aliphatic chain and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"
    
    # Count the number of carboxylic acid groups. A fatty acid should have just one
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
         return False, "Too many carboxylic acid groups"

    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # 2. Check for chain length (4-28 carbons, or more for extended)
    if c_count < 4 :
        return False, "Chain too short to be a fatty acid"

    # 3. Check for aliphatic nature (mostly C and H) and allowed substituents
    allowed_atoms = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]  # H, C, N, O, F, P, S, Cl, Br, I
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Non allowed atoms present"
    
    # 4. Branching Check - allow some branching (<= 3 side chains), avoid cyclic substructures.
    methyl_branch = Chem.MolFromSmarts("[CX4]([H])([H])([H])")
    ethyl_branch = Chem.MolFromSmarts("[CX4]([H])([H])([CX4][H])")
    branch_count = len(mol.GetSubstructMatches(methyl_branch)) + len(mol.GetSubstructMatches(ethyl_branch))
    
    if branch_count > 3:
         return False, "Too many branches for typical fatty acid."

    # Check for rings: do not allow structures with rings other than epoxides
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0 :
        epoxide_pattern = Chem.MolFromSmarts("C1CO1") # allow for epoxide rings
        epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
        if ring_info.NumRings() > len(epoxide_matches):
            for ring in ring_info.AtomRings(): # Check if rings are embedded in the carbon chain, other than epoxide
                is_embedded = True
                for atom_index in ring:
                    atom = mol.GetAtomWithIdx(atom_index)
                    if atom.GetAtomicNum() != 6: # ring should be made of carbons
                         is_embedded = False
                         break;
                    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                    
                    carbon_count = 0
                    for n_idx in neighbors:
                        if mol.GetAtomWithIdx(n_idx).GetAtomicNum() == 6:
                             carbon_count += 1;
                    if carbon_count < 2: # if the carbon is not part of a chain, it's not an issue
                        is_embedded = False
                        break
                if is_embedded:
                    return False, "Cyclic structure (other than epoxide) is not allowed"


    # 5. Check for unsaturation (double or triple bonds)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    has_triple_bond = mol.HasSubstructMatch(triple_bond_pattern)


    # Check for a minimum of 2 carbons between the acid and the unsaturated portion
    # This is because the typical fatty acid has the acid at the end of the chain
    if has_double_bond or has_triple_bond:
        acid_carbon = mol.GetSubstructMatch(carboxylic_acid_pattern)[0]
        unsaturated_match_1 = mol.GetSubstructMatches(double_bond_pattern)
        unsaturated_match_2 = mol.GetSubstructMatches(triple_bond_pattern)

        if unsaturated_match_1:
           for unsat_match in unsaturated_match_1:
                unsaturated_carbon_1 = unsat_match[0]
                unsaturated_carbon_2 = unsat_match[1]
                path_1 = GraphDescriptors.GetShortestPath(mol, acid_carbon, unsaturated_carbon_1) if mol.GetBondBetweenAtoms(acid_carbon,unsaturated_carbon_1) is None else 1
                path_2 = GraphDescriptors.GetShortestPath(mol, acid_carbon, unsaturated_carbon_2) if mol.GetBondBetweenAtoms(acid_carbon,unsaturated_carbon_2) is None else 1
                if path_1 < 2 and path_2 < 2:
                      return False, "Unsaturated bond too close to carboxylic acid"
        if unsaturated_match_2:
            for unsat_match in unsaturated_match_2:
                unsaturated_carbon_1 = unsat_match[0]
                unsaturated_carbon_2 = unsat_match[1]
                path_1 = GraphDescriptors.GetShortestPath(mol, acid_carbon, unsaturated_carbon_1) if mol.GetBondBetweenAtoms(acid_carbon,unsaturated_carbon_1) is None else 1
                path_2 = GraphDescriptors.GetShortestPath(mol, acid_carbon, unsaturated_carbon_2) if mol.GetBondBetweenAtoms(acid_carbon,unsaturated_carbon_2) is None else 1
                if path_1 < 2 and path_2 < 2:
                      return False, "Unsaturated bond too close to carboxylic acid"


    # 6. Molecular weight and complexity
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
       return False, "Molecular weight too low to be a fatty acid"

    if mol_wt > 1000:
       return False, "Molecular weight too high to be a fatty acid"


    return True, "Meets criteria for a fatty acid"