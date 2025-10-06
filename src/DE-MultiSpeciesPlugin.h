#pragma once

#include <Dataset.h>
#include <ViewPlugin.h>
#include <ClusterData/ClusterData.h>
#include <PointData/DimensionPickerAction.h>
#include <PointData/PointData.h>
#include <widgets/DropWidget.h>

#include "AdditionalSettings.h"
#include "ButtonProgressBar.h"
#include "LoadedDatasetsAction.h"
#include "MultiTriggerAction.h"
#include "TableModel.h"
#include "TableSortFilterProxyModel.h"
#include "TableView.h"

#include <array>

#include <QTableWidget>

using namespace mv::plugin;
using namespace mv::gui;
using namespace mv::util;

class QLabel;

class DEMultiSpeciesPlugin : public ViewPlugin
{
    Q_OBJECT

public:

    /**
     * Constructor
     * @param factory Pointer to the plugin factory
     */
    DEMultiSpeciesPlugin(const PluginFactory* factory);

    /** Destructor */
    ~DEMultiSpeciesPlugin() override = default;
    
    /** This function is called by the core after the view plugin has been created */
    void init() override;

    void setPositionDataset(const mv::Dataset<Points>& newPoints);

    void setClustersDataset(const mv::Dataset<Clusters>& newClusters);

public: // Miscellaneous
    /** Get smart pointer to points dataset for point position */
    mv::Dataset<Points>& getPositionDataset() { return _points; }

private:

    // Depending on _clusters set the number of columns in the main table
    void updateTableModel();

    // Calculate min, max and rescale values
    void computeMetaData();

    // if _useSelMapForDE == 2: map IDs from clusters to points via selection map, otherwise use cluster IDs
    const std::vector<unsigned int>& getSpeciesIDs(const size_t species);

public: // Serialization

    /**
    * Load plugin from variant map
    * @param Variant map representation of the plugin
    */
    void fromVariantMap(const QVariantMap& variantMap) override;

    /**
    * Save plugin to variant map
    * @return Variant map representation of the plugin
    */
    QVariantMap toVariantMap() const override;

protected slots:
    void writeToCSV() const;
    void computeDE();
    
    void tableView_clicked(const QModelIndex& index);
    void tableView_selectionChanged(const QItemSelection& selected, const QItemSelection& deselected);

protected:
    using QLabelArray2 = std::array<QLabel, MultiTriggerAction::Size>;

    DropWidget*                             _dropWidget;                /** Widget for drag and drop behavior */
    mv::Dataset<Points>                     _points;                    /** Points smart pointer */
    mv::Dataset<Clusters>                   _clusters;                  /** Clusters smart pointer */
    QLabel*                                 _currentDatasetNameLabel;   /** Label that show the current dataset name */
   
    MultiTriggerAction                      _setSelectionTriggerActions;
    MultiTriggerAction                      _highlightSelectionTriggerActions;
    LoadedDatasetsAction                    _loadedDatasetsAction;
    TriggerAction                           _updateStatisticsAction;
    StringAction                            _filterOnIdAction;
    StringAction                            _selectedIdAction;
    TriggerAction                           _copyToClipboardAction;
    TriggerAction                           _saveToCsvAction;
    TriggerAction                           _openAdditionalSettingsAction;
    DimensionPickerAction                   _currentSelectedDimension;
    AdditionalSettingsDialog                _additionalSettingsDialog;

    QLabelArray2                            _selectedCellsLabel;
    int                                     _totalTableColumns;
    QSharedPointer<TableModel>              _tableItemModel;
    QPointer<TableSortFilterProxyModel>     _sortFilterProxyModel;
    TableView*                              _tableView;
    QPointer<ButtonProgressBar>             _buttonProgressBar;

    QVector<WidgetAction*>                  _serializedActions;
    QByteArray                              _headerState;

    std::vector<QTableWidgetItem*>          _geneTableItems = {};
    std::vector<QTableWidgetItem*>          _diffTableItems = {};

    std::vector<std::vector<float>>         _minValues = {};            // min values for each dimension for each species (cluster)
    std::vector<std::vector<float>>         _rescaleValues = {};        // rescale values for each dimension for each species (cluster)

    std::vector<uint32_t>                   _selectionA = {};
    std::vector<uint32_t>                   _selectionB = {};

    ToggleAction                            _normAction;                // min max normalization
    bool                                    _norm = false;
    int                                     _useSelMapForDE = 0;        // 0: not init, 1: dont, 2: do
    std::vector<unsigned int>               _mappedSpeciesIDs = {};     // if _useSelMapForDE == 2: map IDs from clusters to points via selection map, otherwise use cluster IDs
};


// =============================================================================
// Factory
// =============================================================================
class DEMultiSpeciesPluginFactory : public ViewPluginFactory
{
    Q_INTERFACES(mv::plugin::ViewPluginFactory mv::plugin::PluginFactory)
    Q_OBJECT
    Q_PLUGIN_METADATA(IID   "nl.BioVault.DEMultiSpeciesPlugin"
                      FILE  "PluginInfo.json")

public:

    /** Default constructor */
    DEMultiSpeciesPluginFactory();

    /** Destructor */
    ~DEMultiSpeciesPluginFactory() override {}
    
    /** Creates an instance of the example view plugin */
    ViewPlugin* produce() override;

    /** Returns the data types that are supported by the example view plugin */
    mv::DataTypes supportedDataTypes() const override;

    /**
     * Get plugin trigger actions given \p datasets
     * @param datasets Vector of input datasets
     * @return Vector of plugin trigger actions
     */
    PluginTriggerActions getPluginTriggerActions(const mv::Datasets& datasets) const override;
};
