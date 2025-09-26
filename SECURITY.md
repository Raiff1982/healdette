# Security Policy

## Supported Versions

We maintain security updates for the following versions of Healdette:

| Version | Supported          |
| ------- | ------------------ |
| 2.0.x   | :white_check_mark: |
| 1.0.x   | :x:                |

## Security Features

### Input Validation and Sanitization
- JSON content validation with depth limits
- File type verification using python-magic
- File size restrictions
- Path traversal prevention
- Filename sanitization
- Input format validation

### Access Control
- CSRF protection on all forms
- Rate limiting to prevent abuse
- Secure file operation handling
- Path access restrictions

### Logging and Monitoring
- Security event logging
- Audit logging for configuration changes
- Access attempt tracking
- Validation failure monitoring
- Detailed error logging with sanitization

### Data Protection
- Automated backup system
- File integrity verification
- Secure restoration process
- Version management
- Backup rotation and cleanup

## Reporting a Vulnerability

We take security seriously at Healdette. If you discover a security vulnerability, please follow these steps:

1. **Do Not** disclose the vulnerability publicly
2. Send a detailed report to security@healdette.example.com including:
   - Description of the vulnerability
   - Steps to reproduce
   - Potential impact
   - Suggested fixes (if any)

### What to Expect
- Acknowledgment of your report within 24 hours
- Regular updates on our progress
- Credit for your responsible disclosure (if desired)

### Our Commitments
- Prompt acknowledgment of your report
- Regular updates on progress
- Notification when the vulnerability is fixed
- Credit in our security advisory (unless you prefer to remain anonymous)

## Security Best Practices

### Configuration Security
1. Always validate configuration files before use
2. Use the provided validation tools
3. Keep backups of working configurations
4. Monitor audit logs for unexpected changes

### API Security
1. Use rate limiting for API endpoints
2. Implement proper error handling
3. Validate all inputs
4. Monitor for unusual activity

### File Operation Security
1. Use provided file operation methods
2. Verify file integrity after transfers
3. Maintain regular backups
4. Monitor audit logs

## Security Updates

Security updates are distributed through:
1. GitHub Security Advisories
2. Release Notes
3. Security Notifications

Subscribe to our security notifications by watching this repository.