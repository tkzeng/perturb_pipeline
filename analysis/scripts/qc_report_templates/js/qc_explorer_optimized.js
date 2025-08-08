// QC Report Explorer - Optimized Version with Lazy Loading

document.addEventListener('alpine:init', () => {
    Alpine.data('qcReport', () => ({
        // State
        manifest: null,
        loading: true,
        currentView: 'overview',
        saturationTab: 'gex',
        compareMode: 'metric',
        consolidatedTab: 'sample',
        consolidatedDataSource: 'main_raw',
        collapsedGroups: {},
        
        // Category data (loaded on demand)
        categoryData: {},
        categoryLoading: {},
        
        // Filters
        filters: {
            pool: '',
            sample: '',
            method: '',
            source: '',
            processing: ''
        },
        
        // Computed data
        samples: [],
        generatedDate: '',
        totalSamples: 0,
        plotCounts: {},
        
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
        
        // Consolidated metric groups
        consolidatedMetricGroups: {
            'Cell Counts': ['n_cells', 'n_noncells', 'fraction_cells'],
            'Gene Expression Quality': ['mean_umis_per_cell', 'median_umis_per_cell', 'mean_genes_per_cell', 'median_genes_per_cell', 'mean_umis_per_noncell', 'median_umis_per_noncell', 'mean_genes_per_noncell', 'median_genes_per_noncell'],
            'Guide Metrics': ['guides_per_cell', 'guide_umis_per_cell', 'fraction_cells_with_guides', 'umis_per_guide_per_cell'],
            'Quality Control': ['mean_pct_mito_cells', 'median_pct_mito_cells', 'mean_pct_mito_noncells', 'median_pct_mito_noncells']
        },
        
        // Initialize - now only loads main manifest
        async init() {
            try {
                const response = await fetch('manifest.json');
                this.manifest = await response.json();
                
                // Process basic manifest data
                this.processMainManifest();
                
                // Initialize sample table
                this.$nextTick(() => {
                    this.initializeSampleTable();
                });
                
                // Load category data based on initial view
                await this.loadCategoryData(this.currentView);
                
                this.loading = false;
            } catch (error) {
                console.error('Error loading manifest:', error);
                this.loading = false;
            }
        },
        
        // Process main manifest (lightweight)
        processMainManifest() {
            if (!this.manifest) return;
            
            this.generatedDate = new Date(this.manifest.generated).toLocaleDateString();
            this.samples = this.manifest.samples || [];
            this.totalSamples = this.samples.length;
            this.plotCounts = this.manifest.plot_counts || {};
            
            // Process tables
            this.tables = this.manifest.tables || [];
        },
        
        // Load category-specific data on demand
        async loadCategoryData(category) {
            // Map view names to category names
            const categoryMap = {
                'overview': null,  // No plots for overview
                'cell-calling': 'cell_calling',
                'saturation': 'saturation',
                'per-cell-qc': 'per_cell',
                'compare-samples': 'per_cell',  // Uses same data
                'consolidated': 'consolidated'
            };
            
            const categoryName = categoryMap[category];
            if (!categoryName || this.categoryData[categoryName]) {
                return; // Already loaded or not needed
            }
            
            // Show loading state
            this.categoryLoading[categoryName] = true;
            
            try {
                const response = await fetch(`manifest-${categoryName}.json`);
                const data = await response.json();
                
                // Store the loaded data
                this.categoryData[categoryName] = data.plots || [];
                
                console.log(`Loaded ${data.plot_count} plots for ${categoryName}`);
            } catch (error) {
                console.error(`Error loading ${categoryName} data:`, error);
                this.categoryData[categoryName] = [];
            } finally {
                this.categoryLoading[categoryName] = false;
            }
        },
        
        // Watch for view changes
        async handleViewChange(newView) {
            await this.loadCategoryData(newView);
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
            // Get from loaded category data
            const methods = new Set();
            Object.values(this.categoryData).forEach(plots => {
                plots.forEach(plot => {
                    if (plot.metadata?.method) {
                        methods.add(plot.metadata.method);
                    }
                });
            });
            return [...methods].sort();
        },
        
        // Get plots for current view
        get currentPlots() {
            const categoryMap = {
                'cell-calling': 'cell_calling',
                'saturation': 'saturation',
                'per-cell-qc': 'per_cell',
                'compare-samples': 'per_cell',
                'consolidated': 'consolidated'
            };
            
            const category = categoryMap[this.currentView];
            return this.categoryData[category] || [];
        },
        
        get filteredPlots() {
            return this.currentPlots.filter(plot => {
                const metadata = plot.metadata || {};
                
                if (this.filters.pool && metadata.pool !== this.filters.pool) return false;
                if (this.filters.sample && metadata.sample !== this.filters.sample) return false;
                if (this.filters.method && metadata.method !== this.filters.method) return false;
                if (this.filters.source && metadata.source !== this.filters.source) return false;
                if (this.filters.processing && metadata.processing !== this.filters.processing) return false;
                
                return true;
            });
        },
        
        get plotsGroupedByMetric() {
            const grouped = {};
            this.filteredPlots.forEach(plot => {
                const metric = plot.metadata?.metric;
                if (metric) {
                    if (!grouped[metric]) grouped[metric] = [];
                    grouped[metric].push(plot);
                }
            });
            
            // Sort plots within each metric
            Object.keys(grouped).forEach(metric => {
                grouped[metric].sort((a, b) => {
                    const sampleA = a.metadata?.sample || '';
                    const sampleB = b.metadata?.sample || '';
                    return sampleA.localeCompare(sampleB);
                });
            });
            
            return grouped;
        },
        
        get plotsGroupedByMethod() {
            const grouped = {};
            this.filteredPlots.forEach(plot => {
                const method = plot.metadata?.method;
                if (method) {
                    if (!grouped[method]) grouped[method] = [];
                    grouped[method].push(plot);
                }
            });
            return grouped;
        },
        
        // Methods
        filteredPlotsForType(type) {
            return this.filteredPlots.filter(plot => {
                const category = plot.metadata?.category;
                
                if (type === 'cell_calling') return category === 'cell_calling';
                if (type === 'umi_saturation') return category === 'saturation' && plot.metadata?.metric !== 'guide_umi_saturation';
                if (type === 'guide_umi_saturation') return category === 'saturation' && plot.metadata?.metric === 'guide_umi_saturation';
                if (type === 'consolidated') return category === 'consolidated';
                if (type === 'per_cell') return category === 'per_cell';
                
                return false;
            });
        },
        
        filteredPlotsForMetric(metric) {
            return this.filteredPlots.filter(plot => plot.metadata?.metric === metric);
        },
        
        getMethodsForMetric(metric) {
            const plots = this.plotsGroupedByMetric[metric] || [];
            return [...new Set(plots.map(p => p.metadata?.method).filter(Boolean))].sort();
        },
        
        filterPlotsByMethod(plots, method) {
            if (!method) return plots;
            return plots.filter(p => p.metadata?.method === method);
        },
        
        groupByMetric(plots) {
            const grouped = {};
            plots.forEach(plot => {
                const metric = plot.metadata?.metric || 'unknown';
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
        
        // Get consolidated plots organized by metric group
        getConsolidatedPlotsByGroup(aggregationLevel, dataSource) {
            const consolidated = this.filteredPlotsForType('consolidated');
            
            // Filter by aggregation level and data source
            const filtered = consolidated.filter(plot => {
                const metadata = plot.metadata || {};
                
                const levelMatch = metadata.aggregation_level === aggregationLevel;
                const sourceMatch = (dataSource === 'main_raw' && metadata.source === 'main' && metadata.processing === 'raw') ||
                                  (dataSource === 'all_merged' && metadata.source === 'all' && metadata.processing === 'merged');
                
                return levelMatch && sourceMatch;
            });
            
            // Group by metric type
            const grouped = {};
            filtered.forEach(plot => {
                const metadata = plot.metadata || {};
                const filename = plot.path.split('/').pop().replace('.png', '');
                
                // Extract metric name
                let metricName = filename;
                metricName = metricName.replace(/^(sample|biological_sample|well)_/, '');
                metricName = metricName.replace(/_(?:BarcodeRanks_Inflection|BarcodeRanks_Knee|EmptyDrops_FDR\d+|Expected_Cells|UMI_Threshold)/, '');
                
                // Determine metric group
                let groupName = 'Other Metrics';
                for (const [group, metrics] of Object.entries(this.consolidatedMetricGroups)) {
                    if (metrics.some(m => metricName.includes(m))) {
                        groupName = group;
                        break;
                    }
                }
                
                if (!grouped[groupName]) grouped[groupName] = [];
                
                grouped[groupName].push({
                    ...plot,
                    metricName: metricName,
                    variant: metadata.scale === 'log' ? 'Log Scale' : 'Linear Scale',
                    cellCallingMethod: this.extractCellCallingMethod(filename)
                });
            });
            
            // Sort plots within each group
            Object.keys(grouped).forEach(group => {
                grouped[group].sort((a, b) => {
                    const metricCompare = a.metricName.localeCompare(b.metricName);
                    if (metricCompare !== 0) return metricCompare;
                    
                    const methodCompare = (a.cellCallingMethod || '').localeCompare(b.cellCallingMethod || '');
                    if (methodCompare !== 0) return methodCompare;
                    
                    return a.variant.localeCompare(b.variant);
                });
            });
            
            return grouped;
        },
        
        // Extract cell calling method from filename
        extractCellCallingMethod(filename) {
            const methods = {
                'BarcodeRanks_Inflection': 'Barcode Ranks (Inflection)',
                'BarcodeRanks_Knee': 'Barcode Ranks (Knee)',
                'EmptyDrops_FDR0001': 'EmptyDrops (FDR 0.001)',
                'EmptyDrops_FDR001': 'EmptyDrops (FDR 0.01)',
                'EmptyDrops_FDR005': 'EmptyDrops (FDR 0.05)',
                'Expected_Cells': 'Expected Cells',
                'UMI_Threshold': 'UMI Threshold'
            };
            
            for (const [key, value] of Object.entries(methods)) {
                if (filename.includes(key)) return value;
            }
            return null;
        },
        
        // Group plots into pairs (linear + log variants)
        groupPlotPairs(plots) {
            const pairs = {};
            
            plots.forEach(plot => {
                const baseKey = plot.path.replace('_log.png', '.png');
                
                if (!pairs[baseKey]) pairs[baseKey] = [];
                pairs[baseKey].push(plot);
            });
            
            return Object.values(pairs).map(pair => {
                return pair.sort((a, b) => a.path.includes('_log') ? 1 : -1);
            });
        },
        
        // Toggle metric group collapse state
        toggleMetricGroup(groupName) {
            this.collapsedGroups[groupName] = !this.collapsedGroups[groupName];
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

// Add lazy loading for images using Intersection Observer
document.addEventListener('DOMContentLoaded', () => {
    // Create intersection observer for lazy loading images
    const imageObserver = new IntersectionObserver((entries, observer) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                const img = entry.target;
                if (img.dataset.src && !img.src) {
                    img.src = img.dataset.src;
                    img.onload = () => {
                        img.classList.add('loaded');
                    };
                }
                observer.unobserve(img);
            }
        });
    }, {
        rootMargin: '50px'
    });
    
    // Watch for new images added to DOM
    const observeImages = () => {
        document.querySelectorAll('img[data-src]').forEach(img => {
            if (!img.src) {
                imageObserver.observe(img);
            }
        });
    };
    
    // Initial observation
    setTimeout(observeImages, 100);
    
    // Re-observe when Alpine updates the DOM
    document.addEventListener('alpine:updated', observeImages);
});