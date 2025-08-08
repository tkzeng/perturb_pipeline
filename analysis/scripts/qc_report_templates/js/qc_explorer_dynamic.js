// QC Report Explorer - Fully Dynamic Version
// NO HARDCODING!

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
        expandedSamples: {},
        
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
        
        // Dynamic data from manifest
        samples: [],
        generatedDate: '',
        totalSamples: 0,
        plotCounts: {},
        organization: {},
        
        // Initialize - loads main manifest with organization info
        async init() {
            try {
                const response = await fetch('manifest.json');
                this.manifest = await response.json();
                
                // Process manifest data
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
        
        // Process main manifest with dynamic organization
        processMainManifest() {
            if (!this.manifest) return;
            
            this.generatedDate = new Date(this.manifest.generated).toLocaleDateString();
            this.samples = this.manifest.samples || [];
            this.totalSamples = this.samples.length;
            this.plotCounts = this.manifest.plot_counts || {};
            
            // Dynamic organization from manifest
            this.organization = this.manifest.organization || {};
            this.tables = this.manifest.tables || [];
        },
        
        // Load category-specific data on demand
        async loadCategoryData(category) {
            const categoryMap = {
                'overview': null,
                'cell-calling': 'cell_calling',
                'saturation': 'saturation',
                'per-cell-qc': 'per_cell',
                'compare-samples': 'per_cell',
                'consolidated': 'consolidated'
            };
            
            const categoryName = categoryMap[category];
            if (!categoryName || this.categoryData[categoryName]) {
                return;
            }
            
            this.categoryLoading[categoryName] = true;
            
            try {
                const response = await fetch(`manifest-${categoryName}.json`);
                const data = await response.json();
                
                // Store the loaded data with hierarchy
                this.categoryData[categoryName] = {
                    plots: data.plots || [],
                    hierarchy: data.hierarchy || {}
                };
                
                console.log(`Loaded ${data.plot_count} plots for ${categoryName}`);
            } catch (error) {
                console.error(`Error loading ${categoryName} data:`, error);
                this.categoryData[categoryName] = { plots: [], hierarchy: {} };
            } finally {
                this.categoryLoading[categoryName] = false;
            }
        },
        
        // Watch for view changes
        async handleViewChange(newView) {
            await this.loadCategoryData(newView);
        },
        
        // Computed properties - all dynamic!
        get gexSamples() {
            return this.samples.filter(s => s.sample_type === 'gex');
        },
        
        get guideSamples() {
            return this.samples.filter(s => s.sample_type === 'guide');
        },
        
        get availablePools() {
            return this.organization.pools || [];
        },
        
        get availableSamples() {
            return [...new Set(this.samples.map(s => s.sample_id).filter(Boolean))].sort();
        },
        
        get availableMethods() {
            return this.organization.methods || [];
        },
        
        get availableMetrics() {
            return this.organization.metrics || [];
        },
        
        // Get current category data
        getCategoryData(category) {
            return this.categoryData[category] || { plots: [], hierarchy: {} };
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
            const data = this.getCategoryData(category);
            return data.plots || [];
        },
        
        get currentHierarchy() {
            const categoryMap = {
                'cell-calling': 'cell_calling',
                'saturation': 'saturation',
                'per-cell-qc': 'per_cell',
                'compare-samples': 'per_cell',
                'consolidated': 'consolidated'
            };
            
            const category = categoryMap[this.currentView];
            const data = this.getCategoryData(category);
            return data.hierarchy || {};
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
        
        // Dynamic grouping based on metadata
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
        
        // Get plots organized by sample hierarchy
        getPerCellPlotsBySample() {
            const hierarchy = this.currentHierarchy;
            const filtered = {};
            
            // Apply filters to hierarchy
            Object.entries(hierarchy).forEach(([sample, metrics]) => {
                // Check if sample passes filters
                const sampleMetadata = { sample };
                if (this.filters.sample && !sample.includes(this.filters.sample)) return;
                
                // Filter metrics
                const filteredMetrics = {};
                Object.entries(metrics).forEach(([metric, plots]) => {
                    const filteredPlots = plots.filter(plot => {
                        const metadata = plot.metadata || {};
                        if (this.filters.pool && metadata.pool !== this.filters.pool) return false;
                        if (this.filters.method && metadata.method !== this.filters.method) return false;
                        if (this.filters.source && metadata.source !== this.filters.source) return false;
                        if (this.filters.processing && metadata.processing !== this.filters.processing) return false;
                        return true;
                    });
                    
                    if (filteredPlots.length > 0) {
                        filteredMetrics[metric] = filteredPlots;
                    }
                });
                
                if (Object.keys(filteredMetrics).length > 0) {
                    filtered[sample] = filteredMetrics;
                }
            });
            
            return filtered;
        },
        
        // Get consolidated plots organized dynamically
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
            
            // Group by metric dynamically
            const grouped = {};
            filtered.forEach(plot => {
                const metric = plot.metadata?.metric || 'other';
                const metricGroup = this.getMetricGroup(metric);
                
                if (!grouped[metricGroup]) grouped[metricGroup] = [];
                
                grouped[metricGroup].push({
                    ...plot,
                    metricName: metric,
                    variant: plot.metadata?.scale === 'log' ? 'Log Scale' : 'Linear Scale',
                    cellCallingMethod: plot.metadata?.method
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
        
        // Dynamically determine metric group based on metric name
        getMetricGroup(metric) {
            // Group metrics by patterns in their names
            if (metric.includes('cell') || metric.includes('fraction')) {
                return 'Cell Counts & Fractions';
            } else if (metric.includes('umi') || metric.includes('gene')) {
                return 'Gene Expression Metrics';
            } else if (metric.includes('guide')) {
                return 'Guide Metrics';
            } else if (metric.includes('mito') || metric.includes('quality')) {
                return 'Quality Control Metrics';
            } else {
                return 'Other Metrics';
            }
        },
        
        // Methods
        filteredPlotsForType(type) {
            return this.filteredPlots.filter(plot => {
                const category = plot.metadata?.category;
                return category === type;
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
        
        // Format metric name for display (dynamic)
        formatMetricName(metric) {
            if (!metric) return '';
            
            // Use display name from metadata if available
            const plot = this.currentPlots.find(p => p.metadata?.metric === metric);
            if (plot?.metadata?.metric_display) {
                return plot.metadata.metric_display;
            }
            
            // Otherwise, auto-format
            return metric
                .replace(/_/g, ' ')
                .replace(/\b\w/g, l => l.toUpperCase())
                .replace(/Umis/g, 'UMIs')
                .replace(/Umi /g, 'UMI ')
                .replace(/Pct /g, '% ')
                .replace(/Vs /g, 'vs ');
        },
        
        // Group plots into pairs (linear + log variants)
        groupPlotPairs(plots) {
            const pairs = {};
            
            plots.forEach(plot => {
                // Create a key without the scale info
                const metadata = plot.metadata || {};
                const key = `${metadata.metric}_${metadata.method}`;
                
                if (!pairs[key]) pairs[key] = [];
                pairs[key].push(plot);
            });
            
            // Convert to array and ensure proper pairing
            return Object.values(pairs).map(pair => {
                return pair.sort((a, b) => {
                    // Linear before log
                    const scaleA = a.metadata?.scale || 'linear';
                    const scaleB = b.metadata?.scale || 'linear';
                    return scaleA.localeCompare(scaleB);
                });
            });
        },
        
        // Toggle functions
        toggleMetricGroup(groupName) {
            this.collapsedGroups[groupName] = !this.collapsedGroups[groupName];
        },
        
        toggleSample(sampleId) {
            this.expandedSamples[sampleId] = !this.expandedSamples[sampleId];
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