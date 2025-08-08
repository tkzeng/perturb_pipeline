// QC Report Explorer - Alpine.js Implementation

document.addEventListener('alpine:init', () => {
    Alpine.data('qcReport', () => ({
        // State
        manifest: null,
        loading: true,
        currentView: 'overview',
        saturationTab: 'gex',
        compareMode: 'metric',
        
        // Filters
        filters: {
            pool: '',
            sample: '',
            method: '',
            source: '',
            processing: ''
        },
        
        // Computed data
        plots: [],
        tables: [],
        samples: [],
        generatedDate: '',
        totalSamples: 0,
        
        // Metric definitions
        metricTypes: {
            'umis_per_cell': 'UMIs per Cell',
            'genes_per_cell': 'Genes per Cell',
            'pct_mito': 'Mitochondrial Percentage',
            'pct_counts_mt': 'Mitochondrial Count Percentage',
            'guide_umis_per_cell': 'Guide UMIs per Cell',
            'guides_per_cell': 'Guides per Cell',
            'umi_vs_genes_scatter': 'UMI vs Genes Scatter',
            'total_counts': 'Total Counts',
            'n_genes_by_counts': 'Genes by Counts'
        },
        
        // Initialize
        async init() {
            try {
                const response = await fetch('manifest.json');
                this.manifest = await response.json();
                
                // Process manifest data
                this.processManifest();
                
                // Initialize DataTables after manifest loads
                this.$nextTick(() => {
                    this.initializeSampleTable();
                });
                
                this.loading = false;
            } catch (error) {
                console.error('Error loading manifest:', error);
                this.loading = false;
            }
        },
        
        // Process manifest and extract metadata
        processManifest() {
            if (!this.manifest) return;
            
            this.generatedDate = new Date(this.manifest.generated).toLocaleDateString();
            this.samples = this.manifest.samples || [];
            this.totalSamples = this.samples.length;
            
            // Process plots with metadata extraction
            this.plots = (this.manifest.plots || []).map(plotPath => {
                const metadata = this.extractPlotMetadata(plotPath);
                return {
                    path: plotPath,
                    ...metadata
                };
            });
            
            // Process tables
            this.tables = this.manifest.tables || [];
        },
        
        // Extract metadata from plot filename
        extractPlotMetadata(plotPath) {
            const filename = plotPath.split('/').pop().replace('.png', '');
            const parts = filename.split('_');
            
            // Try to parse different filename patterns
            let metadata = {
                filename: filename,
                title: filename.replace(/_/g, ' ')
            };
            
            // Pattern: pool{n}_{sample_id}_{type}_{metric}_{method}_{source}_{processing}
            // Example: pool1_gex_1_sample_umis_per_cell_EmptyDrops_FDR01_main_raw
            if (plotPath.includes('per_cell/')) {
                // Find the pool and sample
                const poolMatch = filename.match(/pool\d+/);
                const sampleMatch = filename.match(/(gex|guide)_\d+/);
                
                if (poolMatch && sampleMatch) {
                    metadata.pool = poolMatch[0];
                    metadata.sample = sampleMatch[0];
                    
                    // Extract metric - look for known metric patterns
                    for (const metric of Object.keys(this.metricTypes)) {
                        if (filename.includes(metric)) {
                            metadata.metric = metric;
                            break;
                        }
                    }
                    
                    // Extract method - look for known methods
                    const methods = ['EmptyDrops_FDR0001', 'EmptyDrops_FDR001', 'EmptyDrops_FDR005', 
                                   'EmptyDrops_FDR01', 'BarcodeRanks_Knee', 'BarcodeRanks_Inflection', 
                                   'Expected_Cells', 'UMI_Threshold'];
                    for (const method of methods) {
                        if (filename.includes(method)) {
                            metadata.method = method;
                            break;
                        }
                    }
                    
                    // Extract source and processing
                    const sourceMatch = filename.match(/(main|undetermined|all)_(raw|recovered|merged)/);
                    if (sourceMatch) {
                        metadata.source = sourceMatch[1];
                        metadata.processing = sourceMatch[2];
                    }
                }
            }
            // Pattern for cell calling: pool{n}_{sample_id}_barcode_rank_{source}_{processing}
            else if (plotPath.includes('cell_calling/')) {
                const poolMatch = filename.match(/pool\d+/);
                const sampleMatch = filename.match(/(gex|guide)_\d+/);
                
                if (poolMatch && sampleMatch) {
                    metadata.pool = poolMatch[0];
                    metadata.sample = sampleMatch[0];
                    metadata.type = 'cell_calling';
                    
                    const sourceMatch = filename.match(/(main|undetermined|all)_(raw|recovered|merged)/);
                    if (sourceMatch) {
                        metadata.source = sourceMatch[1];
                        metadata.processing = sourceMatch[2];
                    }
                }
            }
            // Pattern for saturation
            else if (plotPath.includes('saturation/')) {
                const poolMatch = filename.match(/pool\d+/);
                const sampleMatch = filename.match(/(gex|guide)_\d+/);
                
                if (poolMatch && sampleMatch) {
                    metadata.pool = poolMatch[0];
                    metadata.sample = sampleMatch[0];
                    metadata.type = filename.includes('guide_umi_saturation') ? 'guide_umi_saturation' : 'umi_saturation';
                }
            }
            // Consolidated plots
            else if (plotPath.includes('consolidated/')) {
                metadata.type = 'consolidated';
                
                const sourceMatch = filename.match(/(main|undetermined|all)_(raw|recovered|merged)/);
                if (sourceMatch) {
                    metadata.source = sourceMatch[1];
                    metadata.processing = sourceMatch[2];
                }
            }
            
            return metadata;
        },
        
        // Computed properties
        get gexSamples() {
            return this.samples.filter(s => s.sample_type === 'gex');
        },
        
        get guideSamples() {
            return this.samples.filter(s => s.sample_type === 'guide');
        },
        
        get availablePools() {
            return [...new Set(this.samples.map(s => s.pool).filter(Boolean))].sort();
        },
        
        get availableSamples() {
            return [...new Set(this.samples.map(s => s.sample_id).filter(Boolean))].sort();
        },
        
        get availableMethods() {
            return [...new Set(this.plots.map(p => p.method).filter(Boolean))].sort();
        },
        
        get filteredPlots() {
            return this.plots.filter(plot => {
                if (this.filters.pool && plot.pool !== this.filters.pool) return false;
                if (this.filters.sample && plot.sample !== this.filters.sample) return false;
                if (this.filters.method && plot.method !== this.filters.method) return false;
                if (this.filters.source && plot.source !== this.filters.source) return false;
                if (this.filters.processing && plot.processing !== this.filters.processing) return false;
                return true;
            });
        },
        
        get plotsGroupedByMetric() {
            const grouped = {};
            this.filteredPlots.forEach(plot => {
                if (plot.metric) {
                    if (!grouped[plot.metric]) grouped[plot.metric] = [];
                    grouped[plot.metric].push(plot);
                }
            });
            // Sort plots within each metric by sample name
            Object.keys(grouped).forEach(metric => {
                grouped[metric].sort((a, b) => (a.sample || '').localeCompare(b.sample || ''));
            });
            return grouped;
        },
        
        get plotsGroupedByMethod() {
            const grouped = {};
            this.filteredPlots.forEach(plot => {
                if (plot.method) {
                    if (!grouped[plot.method]) grouped[plot.method] = [];
                    grouped[plot.method].push(plot);
                }
            });
            return grouped;
        },
        
        // Methods
        filteredPlotsForType(type) {
            return this.filteredPlots.filter(plot => {
                if (type === 'cell_calling') return plot.path.includes('cell_calling/');
                if (type === 'umi_saturation') return plot.path.includes('saturation/') && !plot.path.includes('guide_umi');
                if (type === 'guide_umi_saturation') return plot.path.includes('guide_umi_saturation');
                if (type === 'consolidated') return plot.path.includes('consolidated/');
                return plot.type === type;
            });
        },
        
        filteredPlotsForMetric(metric) {
            return this.filteredPlots.filter(plot => plot.metric === metric);
        },
        
        getMethodsForMetric(metric) {
            const plots = this.plotsGroupedByMetric[metric] || [];
            return [...new Set(plots.map(p => p.method).filter(Boolean))].sort();
        },
        
        filterPlotsByMethod(plots, method) {
            if (!method) return plots;
            return plots.filter(p => p.method === method);
        },
        
        groupByMetric(plots) {
            const grouped = {};
            plots.forEach(plot => {
                const metric = plot.metric || 'unknown';
                if (!grouped[metric]) grouped[metric] = [];
                grouped[metric].push(plot);
            });
            return grouped;
        },
        
        resetFilters() {
            this.filters = {
                pool: '',
                sample: '',
                method: '',
                source: '',
                processing: ''
            };
        },
        
        // Initialize DataTables for sample table
        initializeSampleTable() {
            if (this.currentView === 'overview' && this.samples.length > 0) {
                $('#sample-table').DataTable({
                    data: this.samples,
                    columns: [
                        { data: 'pool', title: 'Pool' },
                        { data: 'sample_id', title: 'Sample ID' },
                        { data: 'sample_type', title: 'Type' },
                        { data: 'expected_cells', title: 'Expected Cells' }
                    ],
                    pageLength: 25,
                    destroy: true
                });
            }
        },
        
        // Load tables dynamically
        loadTable(tablePath, containerId) {
            fetch(tablePath)
                .then(response => response.text())
                .then(data => {
                    const lines = data.trim().split('\n');
                    const headers = lines[0].split('\t');
                    const rows = lines.slice(1).map(line => line.split('\t'));
                    
                    const tableId = 'table-' + Math.random().toString(36).substr(2, 9);
                    const container = document.getElementById(containerId);
                    
                    const tableHtml = `
                        <div class="table-container">
                            <table id="${tableId}" class="display"></table>
                        </div>
                    `;
                    container.innerHTML += tableHtml;
                    
                    // Initialize DataTable
                    $(`#${tableId}`).DataTable({
                        data: rows,
                        columns: headers.map(h => ({ title: h })),
                        pageLength: 25,
                        scrollX: true
                    });
                });
        }
    }));
});