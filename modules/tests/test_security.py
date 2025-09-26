import unittest
from pathlib import Path
import json
import tempfile
import os
from ..security import SecurityValidator
from ..sanitizer import OutputSanitizer
from ..logging_manager import SecurityLogger

class TestSecurityFeatures(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.security = SecurityValidator()
        self.sanitizer = OutputSanitizer()
        self.logger = SecurityLogger(log_dir=self.temp_dir)
        
    def tearDown(self):
        # Clean up temporary files
        for root, dirs, files in os.walk(self.temp_dir):
            for f in files:
                os.unlink(os.path.join(root, f))
        os.rmdir(self.temp_dir)
        
    def test_json_validation(self):
        """Test JSON content validation."""
        # Valid JSON
        valid_json = {
            "name": "test",
            "values": [1, 2, 3],
            "nested": {"key": "value"}
        }
        self.assertTrue(self.security.is_safe_json_content(valid_json))
        
        # Invalid JSON (too deep)
        invalid_json = {}
        current = invalid_json
        for _ in range(25):  # Exceeds MAX_JSON_DEPTH
            current["nested"] = {}
            current = current["nested"]
        with self.assertRaises(ValueError):
            self.security.is_safe_json_content(invalid_json)
            
    def test_filename_validation(self):
        """Test filename validation."""
        self.assertTrue(self.security.is_safe_filename("test.json"))
        self.assertFalse(self.security.is_safe_filename("../test.json"))
        self.assertFalse(self.security.is_safe_filename("test.exe"))
        
    def test_output_sanitization(self):
        """Test output sanitization."""
        # Test HTML sanitization
        html_input = '<script>alert("xss")</script>'
        sanitized = self.sanitizer.sanitize_html(html_input)
        self.assertNotIn("<script>", sanitized)
        
        # Test error message sanitization
        error_msg = "Error in C:\\Users\\test\\file.py: stack trace info"
        sanitized = self.sanitizer.sanitize_error_message(error_msg)
        self.assertNotIn("C:\\Users", sanitized)
        
        # Test JSON sanitization
        json_input = {
            "html": "<img src=x onerror=alert(1)>",
            "nested": {"script": "<script>evil()</script>"}
        }
        sanitized = self.sanitizer.sanitize_json_output(json_input)
        self.assertNotIn("<script>", json.dumps(sanitized))
        
    def test_security_logging(self):
        """Test security logging functionality."""
        # Test security event logging
        self.logger.log_security_event(
            "TEST_EVENT",
            "Test security event",
            {"detail": "test"}
        )
        
        log_file = Path(self.temp_dir) / "security.log"
        self.assertTrue(log_file.exists())
        with open(log_file) as f:
            log_content = f.read()
            self.assertIn("TEST_EVENT", log_content)
            
        # Test audit logging
        self.logger.log_audit_event(
            "TEST_ACTION",
            "test_user",
            "test_resource"
        )
        
        audit_file = Path(self.temp_dir) / "audit.log"
        self.assertTrue(audit_file.exists())
        with open(audit_file) as f:
            audit_content = f.read()
            self.assertIn("TEST_ACTION", audit_content)

if __name__ == '__main__':
    unittest.main()