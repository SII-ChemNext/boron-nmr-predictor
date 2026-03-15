"""Custom exception classes"""


class InvalidSMILESError(Exception):
    """Invalid SMILES exception"""
    pass


class PredictionError(Exception):
    """Exception raised during prediction"""
    pass


class DatabaseError(Exception):
    """Database operation exception"""
    pass


class ValidationError(Exception):
    """Input validation exception"""
    pass
