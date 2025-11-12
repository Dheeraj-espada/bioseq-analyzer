#!/usr/bin/env python3
"""HTML Report Generator"""

import base64
from pathlib import Path
from datetime import datetime
import pandas as pd


class HTMLReportGenerator:
    """Generate HTML reports"""
    
    def __init__(self, project_name="Sequence Analysis"):
        self.project_name = project_name
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.sections = []
        
    def _encode_image(self, image_path):
        """Encode image to base64"""
        try:
            with open(image_path, 'rb') as f:
                encoded = base64.b64encode(f.read()).decode('utf-8')
            return f"data:image/png;base64,{encoded}"
        except:
            return ""
    
    def _generate_css(self):
        """Generate CSS"""
        return """
        <style>
            * { margin: 0; padding: 0; box-sizing: border-box; }
            body {
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.6; color: #333;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                padding: 20px;
            }
            .container {
                max-width: 1200px; margin: 0 auto; background: white;
                border-radius: 10px; box-shadow: 0 10px 40px rgba(0,0,0,0.2);
                overflow: hidden;
            }
            .header {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white; padding: 40px; text-align: center;
            }
            .header h1 { font-size: 2.5em; margin-bottom: 10px; }
            .content { padding: 40px; }
            .section { margin-bottom: 50px; padding-bottom: 30px; border-bottom: 2px solid #eee; }
            .section-title {
                font-size: 2em; color: #667eea; margin-bottom: 20px;
                border-bottom: 3px solid #667eea; padding-bottom: 10px;
            }
            .stats-grid {
                display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 20px; margin: 20px 0;
            }
            .stat-card {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white; padding: 25px; border-radius: 10px;
                box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            }
            .stat-value { font-size: 2.5em; font-weight: bold; margin: 10px 0; }
            .stat-label { font-size: 1em; opacity: 0.9; }
            .visualization { margin: 30px 0; text-align: center; }
            .visualization img {
                max-width: 100%; border-radius: 8px;
                box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            }
            table {
                width: 100%; border-collapse: collapse; margin: 20px 0;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }
            th {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white; padding: 15px; text-align: left;
            }
            td { padding: 12px 15px; border-bottom: 1px solid #eee; }
            tr:hover { background: #f5f5f5; }
            .footer { background: #2c3e50; color: white; text-align: center; padding: 20px; }
            .success-box {
                background: #e8f5e9; border-left: 4px solid #4caf50;
                padding: 15px; margin: 15px 0; border-radius: 4px;
            }
        </style>
        """
    
    def add_summary(self, stats_dict):
        """Add summary section"""
        html = '<div class="section"><h2 class="section-title">📊 Summary</h2><div class="stats-grid">'
        for label, value in stats_dict.items():
            html += f'<div class="stat-card"><div class="stat-label">{label}</div><div class="stat-value">{value}</div></div>'
        html += '</div></div>'
        self.sections.append(html)
    
    def add_table(self, title, dataframe):
        """Add table section"""
        html = f'<div class="section"><h2 class="section-title">{title}</h2>'
        html += dataframe.to_html(index=False, border=0, escape=False)
        html += '</div>'
        self.sections.append(html)
    
    def add_visualizations(self, title, image_paths):
        """Add visualization section"""
        if not isinstance(image_paths, list):
            image_paths = [image_paths]
        
        html = f'<div class="section"><h2 class="section-title">{title}</h2>'
        for img_path in image_paths:
            img_data = self._encode_image(img_path)
            if img_data:
                html += f'<div class="visualization"><img src="{img_data}"></div>'
        html += '</div>'
        self.sections.append(html)
    
    def add_findings(self, findings_list):
        """Add findings section"""
        html = '<div class="section"><h2 class="section-title">🎯 Key Findings</h2>'
        for finding in findings_list:
            html += f'<div class="success-box"><strong>✓</strong> {finding}</div>'
        html += '</div>'
        self.sections.append(html)
    
    def generate(self, output_file):
        """Generate complete HTML report"""
        html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8">
<title>{self.project_name}</title>
{self._generate_css()}
</head><body>
<div class="container">
    <div class="header">
        <h1>🧬 {self.project_name}</h1>
        <p>Comprehensive Sequence Analysis Report</p>
        <p style="font-size: 0.9em; opacity: 0.8;">Generated: {self.timestamp}</p>
    </div>
    <div class="content">
"""
        for section in self.sections:
            html += section
        
        html += f"""
    </div>
    <div class="footer"><p>Generated by BioSeq Analyzer | {self.timestamp}</p></div>
</div></body></html>
"""
        
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html)
            print(f"✅ HTML report: {output_file}")
            return True
        except Exception as e:
            print(f"❌ Error: {e}")
            return False
