import logging
from logging.handlers import RotatingFileHandler
import os
from pathlib import Path
import json
from datetime import datetime
from typing import Any, Dict, Optional

class SecurityLogger:
    def __init__(self, log_dir: str = "logs"):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(exist_ok=True)
        
        # Set up security event logger
        self.security_logger = logging.getLogger('security')
        self.security_logger.setLevel(logging.INFO)
        
        # Create rotating file handler for security events
        security_handler = RotatingFileHandler(
            self.log_dir / 'security.log',
            maxBytes=1024*1024,  # 1MB
            backupCount=10
        )
        security_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        self.security_logger.addHandler(security_handler)
        
        # Set up audit logger
        self.audit_logger = logging.getLogger('audit')
        self.audit_logger.setLevel(logging.INFO)
        
        # Create rotating file handler for audit events
        audit_handler = RotatingFileHandler(
            self.log_dir / 'audit.log',
            maxBytes=1024*1024,  # 1MB
            backupCount=10
        )
        audit_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        self.audit_logger.addHandler(audit_handler)
        
    def log_security_event(self, 
                          event_type: str, 
                          message: str, 
                          details: Optional[Dict[str, Any]] = None,
                          level: str = 'INFO'):
        """Log a security-related event."""
        log_entry = {
            'timestamp': datetime.utcnow().isoformat(),
            'event_type': event_type,
            'message': message,
            'details': details or {}
        }
        
        log_method = getattr(self.security_logger, level.lower())
        log_method(json.dumps(log_entry))
        
    def log_audit_event(self, 
                       action: str, 
                       user: str, 
                       resource: str, 
                       details: Optional[Dict[str, Any]] = None,
                       status: str = 'SUCCESS'):
        """Log an audit event for configuration changes."""
        log_entry = {
            'timestamp': datetime.utcnow().isoformat(),
            'action': action,
            'user': user,
            'resource': resource,
            'status': status,
            'details': details or {}
        }
        
        self.audit_logger.info(json.dumps(log_entry))
        
    def log_validation_failure(self, 
                             validation_type: str, 
                             input_data: Any, 
                             reason: str):
        """Log validation failures."""
        self.log_security_event(
            event_type='VALIDATION_FAILURE',
            message=f'Validation failed for type: {validation_type}',
            details={
                'reason': reason,
                'input': str(input_data)[:200]  # Truncate long inputs
            },
            level='WARNING'
        )
        
    def log_access_attempt(self, 
                          resource: str, 
                          access_type: str, 
                          status: str,
                          client_info: Optional[Dict[str, str]] = None):
        """Log access attempts to resources."""
        self.log_security_event(
            event_type='ACCESS_ATTEMPT',
            message=f'{access_type} access attempt to {resource}',
            details={
                'status': status,
                'client_info': client_info or {}
            },
            level='INFO' if status == 'SUCCESS' else 'WARNING'
        )