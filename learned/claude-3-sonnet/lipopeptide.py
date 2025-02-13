"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import numpy as np
from sklearn.ensemble import RandomForestClassifier

def calculate_descriptors(mol):
    """
    Calculate relevant molecular descriptors for a given molecule.

    Args:
        mol (Mol): RDKit molecule object

    Returns:
        list: List of descriptor values
    """
    descriptors = []
    
    # Calculate molecular weight
    descriptors.append(rdMolDescriptors.CalcExactMolWt(mol))
    
    # Calculate number of atoms and bonds
    descriptors.append(mol.GetNumAtoms())
    descriptors.append(mol.GetNumBonds())
    
    # Calculate number of rotatable bonds
    descriptors.append(rdMolDescriptors.CalcNumRotatableBonds(mol))
    
    # Calculate topological polar surface area
    descriptors.append(rdMolDescriptors.CalcTPSA(mol))
    
    # Calculate Wiener index
    descriptors.append(rdMolDescriptors.CalcWienerIndex(mol))
    
    return descriptors

def is_lipopeptide(smiles):
    """
    Determine if a molecule is a lipopeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - Classification result and reason
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular descriptors
    descriptors = calculate_descriptors(mol)
    
    # Predict using pre-trained random forest model
    prediction = model.predict([descriptors])
    
    if prediction[0]:
        return True, "Molecule classified as lipopeptide by the model"
    else:
        return False, "Molecule not classified as lipopeptide by the model"

# Load or train the random forest model
X_train, y_train = load_training_data()  # Replace with code to load training data
model = RandomForestClassifier()
model.fit(X_train, y_train)