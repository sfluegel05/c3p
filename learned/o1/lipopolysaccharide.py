"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: lipopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide consists of a lipid A moiety (disaccharide of glucosamine with fatty acids),
    a core polysaccharide (including KDO and heptoses), and an O-antigen polysaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize feature counters
    features_found = []
    criteria_met = 0
    essential_criteria_met = 0

    # 1. Check for lipid A moiety (glucosamine disaccharide with fatty acids)
    # Simplifying, check for at least two glucosamine units connected and substituted with long-chain fatty acids

    # SMARTS for glucosamine unit
    glucosamine_smarts = '[C@H]1([O])[C@@H]([O])[C@H](N)[C@@H]([O])[C@H]([O])[C@@H]1[O]'
    glucosamine = Chem.MolFromSmarts(glucosamine_smarts)
    glucosamine_matches = mol.GetSubstructMatches(glucosamine)
    num_glucosamines = len(glucosamine_matches)

    if num_glucosamines >= 2:
        # Check if connected
        # Create a subgraph of the glucosamine units and check for connectivity
        glucosamine_atoms = [match[0] for match in glucosamine_matches]
        are_connected = False
        for i in range(len(glucosamine_atoms)):
            for j in range(i+1, len(glucosamine_atoms)):
                path = Chem.rdmolops.GetShortestPath(mol, glucosamine_atoms[i], glucosamine_atoms[j])
                if path and len(path) <= 10:  # Arbitrary path length limit
                    are_connected = True
                    break
            if are_connected:
                break
        if are_connected:
            # Now check for fatty acid substitutions on glucosamine units
            fatty_acid_smarts = '[NX3][$([C](=O)[C][C][C][C][C][C][C][C][C][C][C])]'  # Amide-linked fatty acid C12 or longer
            ester_fatty_acid_smarts = '[OX2H][$([C](=O)[C][C][C][C][C][C][C][C][C][C][C])]'  # Ester-linked fatty acid
            fatty_acid = Chem.MolFromSmarts(fatty_acid_smarts)
            ester_fatty_acid = Chem.MolFromSmarts(ester_fatty_acid_smarts)
            num_fatty_acids = len(mol.GetSubstructMatches(fatty_acid)) + len(mol.GetSubstructMatches(ester_fatty_acid))
            if num_fatty_acids >= 4:
                essential_criteria_met +=1
                criteria_met +=1
                features_found.append(f"Contains lipid A moiety with {num_fatty_acids} fatty acid chains")
            else:
                features_found.append(f"Glucosamine disaccharide found but insufficient fatty acid chains ({num_fatty_acids})")
        else:
            features_found.append("Glucosamine units found but not connected as in lipid A")
    else:
        features_found.append(f"Contains {num_glucosamines} glucosamine units, less than required for lipid A")

    # 2. Check for KDO units (octulosonic acid)
    kdo_smarts = '[C@H]1([O])[C@H](C(=O)[O])[C@H]([O])[C@@H]([O])[C@H]([O])[C@@H]1[O]'  # Simplified KDO ring
    kdo = Chem.MolFromSmarts(kdo_smarts)
    kdo_matches = mol.GetSubstructMatches(kdo)
    num_kdo = len(kdo_matches)

    if num_kdo >= 1:
        essential_criteria_met +=1
        criteria_met +=1
        features_found.append(f"Contains {num_kdo} KDO (octulosonic acid) units")
    else:
        features_found.append("No KDO (octulosonic acid) units found")

    # 3. Check for heptose units
    heptose_smarts = '[C@H]1([O])[C@H]([O])[C@H]([O])[C@H]([O])[C@H]([O])[C@H]([O])[C@H]1[O]'  # Heptose ring
    heptose = Chem.MolFromSmarts(heptose_smarts)
    heptose_matches = mol.GetSubstructMatches(heptose)
    num_heptose = len(heptose_matches)

    if num_heptose >= 2:
        essential_criteria_met +=1
        criteria_met +=1
        features_found.append(f"Contains {num_heptose} heptose units")
    else:
        features_found.append(f"Contains {num_heptose} heptose units, less than required")

    # 4. Check for multiple sugar units (total sugar rings)
    pyranose_smarts = '[C@H]1([O])[C@@H]([O])[C@H]([O])[C@@H]([O])[C@H]([O])[C@@H]1[O]'  # General pyranose ring
    pyranose = Chem.MolFromSmarts(pyranose_smarts)
    sugar_matches = mol.GetSubstructMatches(pyranose)
    num_sugars = len(sugar_matches)

    if num_sugars >= 5:
        criteria_met +=1
        features_found.append(f"Contains {num_sugars} sugar units")
    else:
        features_found.append(f"Contains {num_sugars} sugar units, less than required")

    # 5. Check Molecular Weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt >= 1000:
        criteria_met +=1
        features_found.append(f"Molecular weight is {mol_wt:.2f} Da")
    else:
        features_found.append(f"Molecular weight is {mol_wt:.2f} Da, less than required")

    # Determine if criteria are met
    if essential_criteria_met >= 2:
        return True, "; ".join(features_found)
    else:
        return False, "; ".join(features_found)