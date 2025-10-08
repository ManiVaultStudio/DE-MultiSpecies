#include "DE-MultiSpeciesPlugin.h"

#include <DatasetsMimeData.h>
#include <util/LearningCenterTutorial.h>

#include <QDebug>
#include <QFile>
#include <QFileDialog>
#include <QJsonArray>
#include <QMimeData>
#include <QPushButton>

#include "AdditionalSettings.h"
#include "WordWrapHeaderView.h"
#include "TutorialUtils.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <iterator>
#include <numeric>
#include <vector>

Q_PLUGIN_METADATA(IID "nl.BioVault.DEMultiSpeciesPlugin")

using namespace mv;

namespace local
{
    static bool is_valid_QByteArray(const QByteArray& state)
    {
        QByteArray data = state;
        QDataStream stream(&data, QIODevice::ReadOnly);
        const int dataStreamVersion = QDataStream::Qt_5_0;
        stream.setVersion(dataStreamVersion);
        int marker;
        int ver;
        stream >> marker;
        stream >> ver;
        bool result = (stream.status() == QDataStream::Ok && (ver == 0));
        return result;
    }

    template <typename T>
    float fround(T n, int d)
    {
        assert(!std::numeric_limits<T>::is_integer); // this function should not be called on integer types
        return static_cast<float>(floor(n * pow(10., d) + 0.5) / pow(10., d));
    }
    template <typename T>
    void resizeNestedVec(std::vector<std::vector<T>>& vecOfVecs, size_t outterSize, size_t innerSize, T val) {
        vecOfVecs.resize(outterSize);
        for (std::vector<T>& innerVec : vecOfVecs) {
            innerVec.resize(innerSize, val);
        }
    }
    template <typename T>
    void sortAndUnique(std::vector<T>& selection) {
        std::sort(selection.begin(), selection.end());
        const auto last = std::unique(selection.begin(), selection.end());
        selection.erase(last, selection.end());
        };
    template<typename T>
    bool is_exact_type(const QVariant& variant)
    {
        auto variantType = variant.metaType();
        auto requestedType = QMetaType::fromType<T>();
        return (variantType == requestedType);
    }
    template<typename T>
    T get_strict_value(const QVariant& variant)
    {
        if (is_exact_type<T>(variant))
            return variant.value<T>();
        else
        {
#ifdef _DEBUG
            qDebug() << "ClusterDEMultiSpeciesPlugin: Error: requested " << QMetaType::fromType<T>().name() << " but value is of type " << variant.metaType().name();
#endif
            return T();
        }
    }
    template <typename FunctionObject>
    void visitAllElements(Dataset<Points>& points, FunctionObject functionObject)
    {
        const auto numDimensions = points->getNumDimensions();
        const auto numRows = points->getNumPoints();
        points->visitData([numRows, numDimensions, functionObject](auto data)
            {
                for (std::size_t r = 0; r < numRows; ++r)
                {
                    for (std::size_t column = 0; column < numDimensions; ++column)
                    {
                        const auto value = data[r][column];
                        functionObject(r, column, value);
                    }
                }
            });
    }
    template <typename RowRange, typename FunctionObject>
    void visitElements(Dataset<Points>& points, const RowRange& rows, FunctionObject functionObject)
    {
        const auto numDimensions = points->getNumDimensions();
        points->visitData([&rows, numDimensions, functionObject](auto data)
            {
                const size_t numRows = rows.size();
                for (std::size_t localRow = 0; localRow < numRows; ++localRow)
                {
                    const auto globalRow = rows[localRow];
                    for (std::size_t column = 0; column < numDimensions; ++column)
                    {
                        const auto value = data[globalRow][column];
                        functionObject(globalRow, localRow, column, value);
                    }
                }
            });
    }

}

DEMultiSpeciesPlugin::DEMultiSpeciesPlugin(const PluginFactory* factory) :
    ViewPlugin(factory),
    _loadedDatasetsAction(this, "Current dataset"),
    _dropWidget(nullptr),
    _points(),
    _clusters(),
    _clustersDataGUID(this, "clusterDataGUID"),
    _currentDatasetNameLabel(new QLabel()),
    _filterOnIdAction(this, "Filter on Id"),
    _selectedIdAction(this, "Last selected Id"),
    _updateStatisticsAction(this, "Calculate Differential Expression"),
    _setSelectionTriggerActions(this, "Set selection triggers", "Set selection %1"),
    _highlightSelectionTriggerActions(this, "Highlight selection triggers", "Highlight selection %1"),
    _sortFilterProxyModel(new TableSortFilterProxyModel),
    _totalTableColumns(0),
    _tableItemModel(new TableModel(nullptr, false)),
    _tableView(nullptr),
    _buttonProgressBar(nullptr),
    _copyToClipboardAction(&getWidget(), "Copy"),
    _saveToCsvAction(&getWidget(), "Save As..."),
    _openAdditionalSettingsAction(&getWidget(), "Open additional settings"),
    _normAction(&getWidget(), "Min-max normalization"),
    _currentSelectedDimension(this, "Selected dimension"),
    _additionalSettingsDialog()
{
    // This line is mandatory if drag and drop behavior is required
    _currentDatasetNameLabel->setAcceptDrops(true);

    // Align text in the center
    _currentDatasetNameLabel->setAlignment(Qt::AlignCenter);

    _normAction.setToolTip("Rescale the data: (selection_mean - global_min) / (global_max - global_min)");

    { // save to CSV

        //addTitleBarMenuAction(&_saveToCsvAction);
        _saveToCsvAction.setIcon(mv::util::StyledIcon("file-csv"));
        _saveToCsvAction.setShortcut(tr("Ctrl+S"));
        _saveToCsvAction.setShortcutContext(Qt::WidgetWithChildrenShortcut);

        connect(&_saveToCsvAction, &TriggerAction::triggered, this, [this]() -> void {
            writeToCSV();
            });
    }

    { // copy to Clipboard

       // addTitleBarMenuAction(&_copyToClipboardAction);
        _copyToClipboardAction.setIcon(mv::util::StyledIcon("copy"));
        _copyToClipboardAction.setShortcut(tr("Ctrl+C"));
        _copyToClipboardAction.setShortcutContext(Qt::WidgetWithChildrenShortcut);

        connect(&_copyToClipboardAction, &TriggerAction::triggered, this, [this]() -> void {
            _tableItemModel->copyToClipboard();
            });
    }

    { // additional settings dialog

        _openAdditionalSettingsAction.setIcon(mv::util::StyledIcon("gears"));
        _additionalSettingsDialog.setCurrentData(_points);

        connect(&_openAdditionalSettingsAction, &TriggerAction::triggered, this, [this]() -> void {
            _additionalSettingsDialog.setCurrentData(_points);
            _additionalSettingsDialog.show();
            });
    }

    _sortFilterProxyModel->setSourceModel(_tableItemModel.get());
    _filterOnIdAction.setSearchMode(true);
    _filterOnIdAction.setClearable(true);
    _filterOnIdAction.setPlaceHolderString("Filter by ID");

    connect(&_updateStatisticsAction, &TriggerAction::triggered, [this](const bool& var)
        {
            _tableItemModel->invalidate();
        });

    connect(&_filterOnIdAction, &mv::gui::StringAction::stringChanged, _sortFilterProxyModel, &TableSortFilterProxyModel::nameFilterChanged);

    connect(&_updateStatisticsAction, &mv::gui::TriggerAction::triggered, this, &DEMultiSpeciesPlugin::computeDE);

    connect(&_normAction, &mv::gui::ToggleAction::toggled, this, [this]()
        {
            if (_normAction.isChecked())
                _norm = true;
            else
                _norm = false;

            _updateStatisticsAction.trigger();
        });

    _serializedActions.append(&_loadedDatasetsAction);
    _serializedActions.append(&_selectedIdAction);
    _serializedActions.append(&_filterOnIdAction);
    _serializedActions.append(&_copyToClipboardAction);
    _serializedActions.append(&_saveToCsvAction);
    _serializedActions.append(&_updateStatisticsAction);
    _serializedActions.append(&_setSelectionTriggerActions);
    _serializedActions.append(&_highlightSelectionTriggerActions);
    _serializedActions.append(&_currentSelectedDimension);
    _serializedActions.append(&_openAdditionalSettingsAction);
}

void DEMultiSpeciesPlugin::init()
{
    QWidget& mainWidget = getWidget();
    _loadedDatasetsAction.initialize(this);

    // Create layout
    auto layout = new QVBoxLayout();

    layout->setContentsMargins(0, 0, 0, 0);

    layout->addWidget(_currentDatasetNameLabel);

    { // toolbar
        QWidget* filterWidget = _filterOnIdAction.createWidget(&mainWidget);
        filterWidget->setContentsMargins(0, 3, 0, 3);

        QHBoxLayout* toolBarLayout = new QHBoxLayout;
        toolBarLayout->addWidget(filterWidget, 8);
        toolBarLayout->addWidget(_normAction.createWidget(&mainWidget), 2);

        layout->addLayout(toolBarLayout);
    }

    { // table view

        _tableView = new TableView(&mainWidget);
        _tableView->setModel(_sortFilterProxyModel);
        _tableView->setSortingEnabled(true);
        _tableView->setSelectionMode(QAbstractItemView::SingleSelection);
        _tableView->setSelectionBehavior(QAbstractItemView::SelectRows);
        _tableView->setContextMenuPolicy(Qt::ActionsContextMenu);

        connect(_tableView->selectionModel(), &QItemSelectionModel::selectionChanged, this, &DEMultiSpeciesPlugin::tableView_selectionChanged);

        WordWrapHeaderView* horizontalHeader = new WordWrapHeaderView(Qt::Horizontal, _tableView, true);

        horizontalHeader->setFirstSectionMovable(false);
        horizontalHeader->setSectionsMovable(true);
        horizontalHeader->setSectionsClickable(true);
        horizontalHeader->sectionResizeMode(QHeaderView::Stretch);
        horizontalHeader->setSectionResizeMode(QHeaderView::Stretch);
        horizontalHeader->setStretchLastSection(true);
        horizontalHeader->setSortIndicator(1, Qt::AscendingOrder);
        horizontalHeader->setDefaultAlignment(Qt::AlignBottom | Qt::AlignLeft | Qt::Alignment(Qt::TextWordWrap));
        _tableView->setHorizontalHeader(horizontalHeader);
        layout->addWidget(_tableView);

        _tableView->addAction(&_saveToCsvAction);
        _tableView->addAction(&_copyToClipboardAction);
        _tableView->addAction(&_openAdditionalSettingsAction);
    }

    {// Progress bar and update button

        _buttonProgressBar = new ButtonProgressBar(&mainWidget, _updateStatisticsAction.createWidget(&mainWidget));
        _buttonProgressBar->setContentsMargins(0, 3, 0, 3);
        _buttonProgressBar->setProgressBarText("No Data Available");
        _buttonProgressBar->setButtonText("Calculate Differential Expression", Qt::black);

        connect(_tableItemModel.get(), &TableModel::statusChanged, _buttonProgressBar, &ButtonProgressBar::showStatus);

        layout->addWidget(_buttonProgressBar);
    }

    updateTableModel();

    // Apply the layout
    getWidget().setLayout(layout);

    // Instantiate new drop widget
    _dropWidget = new DropWidget(_currentDatasetNameLabel);

    // Set the drop indicator widget (the widget that indicates that the view is eligible for data dropping)
    _dropWidget->setDropIndicatorWidget(new DropWidget::DropIndicatorWidget(&getWidget(), "Drag an item points and clusters", ""));

    // Initialize the drop regions
    _dropWidget->initialize([this](const QMimeData* mimeData) -> DropWidget::DropRegions {
        // A drop widget can contain zero or more drop regions
        DropWidget::DropRegions dropRegions;

        const auto datasetsMimeData = dynamic_cast<const DatasetsMimeData*>(mimeData);

        if (datasetsMimeData == nullptr)
            return dropRegions;

        if (datasetsMimeData->getDatasets().count() > 1)
            return dropRegions;

        // Gather information to generate appropriate drop regions
        const auto dataset = datasetsMimeData->getDatasets().first();
        const auto datasetGuiName = dataset->getGuiName();
        const auto datasetId = dataset->getId();
        const auto dataType = dataset->getDataType();
        const auto dataTypes = DataTypes({ PointType, ClusterType });

        // Visually indicate if the dataset is of the wrong data type and thus cannot be dropped
        if (!dataTypes.contains(dataType)) {
            dropRegions << new DropWidget::DropRegion(this, "Incompatible data", "This type of data is not supported", "exclamation-circle", false);
            qDebug() << "ClusterDEMultiSpeciesPlugin: Incompatible data: This type of data is not supported";
        }
        else
        {
            // Accept points datasets drag and drop
            if (dataType == PointType) {

                const auto candidateDataset = mv::data().getDataset<Points>(datasetId);

                if (_points == candidateDataset) { // Dataset cannot be dropped because it is already loaded
                    dropRegions << new DropWidget::DropRegion(this, "Warning", "Data already loaded", "exclamation-circle", false);
                    qDebug() << "ClusterDEMultiSpeciesPlugin: Warning: Data already loaded";
                }
                else { // Dataset can be dropped
                    dropRegions << new DropWidget::DropRegion(this, "Points", QString("Load %1 as points").arg(datasetGuiName), "map-marker-alt", true, [this, candidateDataset]() {
                        setPositionDataset(candidateDataset);
                        });
                }
            }
            else if (dataType == ClusterType) {

                const auto candidateDataset = mv::data().getDataset<Clusters>(datasetId);

                if (_clusters == candidateDataset) { // Dataset cannot be dropped because it is already loaded
                    dropRegions << new DropWidget::DropRegion(this, "Warning", "Data already loaded", "exclamation-circle", false);
                    qDebug() << "ClusterDEMultiSpeciesPlugin: Warning: Data already loaded";
                }
                else { // Dataset can be dropped
                    dropRegions << new DropWidget::DropRegion(this, "Clusters", QString("Load %1 as clusters").arg(datasetGuiName), "map-marker-alt", true, [this, candidateDataset]() {
                        setClustersDataset(candidateDataset);
                        });
                }

            }
        }

        return dropRegions;
        });

    // Respond when the name of the dataset in the dataset reference changes
    connect(&_points, &Dataset<Points>::guiNameChanged, this, [this]() {

        auto newDatasetName = _points->getGuiName();

        // Update the current dataset name label
        _currentDatasetNameLabel->setText(QString("Current points dataset: %1").arg(newDatasetName));

        // Only show the drop indicator when nothing is loaded in the dataset reference
        _dropWidget->setShowDropIndicator(newDatasetName.isEmpty());
        });

    for (std::size_t i = 0; i < _selectedCellsLabel.size(); ++i)
    {
        _selectedCellsLabel[i].setText(QString("(%1 items)").arg(0));
        _selectedCellsLabel[i].setAlignment(Qt::AlignHCenter);
    }

    auto updateSelectionIndices = [this](std::vector<uint32_t>& selection, const QString& selectionName, QLabel& label) {
        if (!_points.isValid())
            return;

        selection = _points->getSelectionIndices();

        // ensure selection is unique and sorted
        local::sortAndUnique(selection);

        label.setText(QString("(%1 items)").arg(selection.size()));

        const auto otherData     = _additionalSettingsDialog.getSelectionMappingSourcePicker().getCurrentDataset<Points>();
        auto& otherDataSelection = _additionalSettingsDialog.getSelection(selectionName);
        otherDataSelection       = otherData.isValid() ? otherData->getSelection<Points>()->indices : std::vector<uint32_t>{};

        local::sortAndUnique(otherDataSelection);

        qDebug() << "ClusterDEMultiSpeciesPlugin: Saved selection " << selectionName << " with " << selection.size() << " items.";

        if (_selectionA.size() != 0 && _selectionB.size() != 0)
            _buttonProgressBar->showStatus(TableModel::Status::OutDated);

        };

    auto highlightSelectionIndices = [this](const std::vector<uint32_t>& selection, const QString& selectionName) {
        if (!_points.isValid())
            return;

        auto otherData = _additionalSettingsDialog.getSelectionMappingSourcePicker().getCurrentDataset<Points>();

        if (otherData.isValid()) {

            // Check if the selection mapping makes sense
            const auto [selectionMapping, numPointsTarget] = getSelectionMappingOtherToCurrent(otherData, _points);
            const bool useOtherSelection = isSurjectiveMappingValid(selectionMapping, numPointsTarget, _points);

            otherData->getSelection<Points>()->indices = useOtherSelection ? _additionalSettingsDialog.getSelection(selectionName) : std::vector<uint32_t>{};
            
            events().notifyDatasetDataSelectionChanged(otherData);
        }
        else {
            _points->setSelectionIndices(selection);
            events().notifyDatasetDataSelectionChanged(_points);
        }

        };

    connect(_setSelectionTriggerActions.getTriggerAction(0), &TriggerAction::triggered, [this, updateSelectionIndices](){
            updateSelectionIndices(_selectionA, "A", _selectedCellsLabel[0]);
        });

    connect(_setSelectionTriggerActions.getTriggerAction(1), &TriggerAction::triggered, [this, updateSelectionIndices](){
            updateSelectionIndices(_selectionB, "B", _selectedCellsLabel[1]);
        });

    connect(_highlightSelectionTriggerActions.getTriggerAction(0), &TriggerAction::triggered, [this, highlightSelectionIndices](){
        highlightSelectionIndices(_selectionA, "A");
        });

    connect(_highlightSelectionTriggerActions.getTriggerAction(1), &TriggerAction::triggered, [this, highlightSelectionIndices](){
            highlightSelectionIndices(_selectionB, "B");
        });

    QGridLayout* selectionLayout = new QGridLayout();
    for (uint8_t i = 0; i < _selectedCellsLabel.size(); ++i)
    {
        selectionLayout->addWidget(_setSelectionTriggerActions.getTriggerAction(i)->createWidget(&mainWidget), 0, i);
        selectionLayout->addWidget(_highlightSelectionTriggerActions.getTriggerAction(i)->createWidget(&mainWidget), 1, i);
        selectionLayout->addWidget(&_selectedCellsLabel[i], 2, i);
    }

    layout->addLayout(selectionLayout);
}

void DEMultiSpeciesPlugin::updateTableModel() {

    _tableItemModel->invalidate();

    int currentColumn = 0;

    if (!_clusters.isValid()) {
        _totalTableColumns = 4;

        _tableItemModel->startModelBuilding(_totalTableColumns, 0);
        _tableItemModel->setHorizontalHeader(currentColumn++, QString("ID"));
        _tableItemModel->setHorizontalHeader(currentColumn++, QString("DE"));
        _tableItemModel->setHorizontalHeader(currentColumn++, QString("Mean (Sel. 1)"));
        _tableItemModel->setHorizontalHeader(currentColumn++, QString("Mean (Sel. 2)"));
        _tableItemModel->endModelBuilding();

        return;
    }

    const auto& clusterNames = _clusters->getClusterNames();
    const auto numSpecies    = clusterNames.size();
       
    _totalTableColumns  = static_cast<int>(1 + numSpecies * 3);

    _tableItemModel->startModelBuilding(_totalTableColumns, 0);
    _tableItemModel->setHorizontalHeader(currentColumn++, QString("ID"));

    for (const auto& speciesName : clusterNames) {
        const auto shortName = speciesName.first(3);
        _tableItemModel->setHorizontalHeader(currentColumn++, QString("DE (%1)").arg(speciesName));
        _tableItemModel->setHorizontalHeader(currentColumn++, QString("Mean 1 (%1)").arg(shortName));
        _tableItemModel->setHorizontalHeader(currentColumn++, QString("Mean 2 (%1)").arg(shortName));
    }

    _tableItemModel->endModelBuilding();
}

void DEMultiSpeciesPlugin::setPositionDataset(const mv::Dataset<Points>& newPoints, bool resetClusterGUID)
{
    if (!newPoints.isValid())
    {
        qDebug() << "DEMultiSpeciesPlugin Warning: invalid points dataset!";
        return;
    }

    _points = newPoints;
    _clusters = {};
    _useSelMapForDE = 0;

    if (resetClusterGUID) {
        _clustersDataGUID.setString("");
    }

    // Update the current data model and dimension picker
    updateTableModel();
    _currentSelectedDimension.setPointsDataset(_points);

    // Only show the drop indicator when nothing is loaded in the dataset reference
    if (_points.isValid() && _clusters.isValid()) {
        _currentDatasetNameLabel->setText(QString("Current datasets: %1 (%2)").arg(_points->getGuiName(), _clusters->getGuiName()));
        _dropWidget->setShowDropIndicator(false);
        computeMetaData();
    }
}

void DEMultiSpeciesPlugin::setClustersDataset(const mv::Dataset<Clusters>& newClusters)
{
    if (!newClusters.isValid())
    {
        qDebug() << "DEMultiSpeciesPlugin Warning: invalid cluster dataset!";
        return;
    }

    if (!_points.isValid())
    {
        qDebug() << "DEMultiSpeciesPlugin Warning: set points dataset first!";
        return;
    }

    // We accept these cases: 
    // 1. the number of points in clusters equals number of points
    // 2. there is a selection map between

    const size_t numPoints      = _points->getNumPoints();
    const size_t numClusterIDs  = std::accumulate(
        newClusters->getClusters().begin(), newClusters->getClusters().end(),
        0ULL,
        [](auto sum, const auto& v) { return sum + v.getIndices().size(); }
    );

    if (numPoints == numClusterIDs) {
        _useSelMapForDE = 1;
        qDebug() << "DEMultiSpeciesPlugin: Points and clusters data cover same number of IDs";
    }
    else if (const auto otherData = _additionalSettingsDialog.getSelectionMappingSourcePicker().getCurrentDataset<Points>(); otherData.isValid()) {
        _useSelMapForDE = 2;
        qDebug() << "DEMultiSpeciesPlugin: Use selection mapping from " << otherData->getGuiName() << " to connect points and clusters";
    }
    else {
        qDebug() << "DEMultiSpeciesPlugin Warning: Points and clusters do not match. Maybe pick a selection mapping data set first and then add the cluster data again.";
        _useSelMapForDE = 0;
        return;
    }

    _clusters = newClusters;
    _clustersDataGUID.setString(_clusters->getId());

    // Update the current data model
    updateTableModel();

    // Only show the drop indicator when nothing is loaded in the dataset reference
    if (_points.isValid() && _clusters.isValid()) {
        _currentDatasetNameLabel->setText(QString("Current datasets: %1 (%2)").arg(_points->getGuiName(), _clusters->getGuiName()));
        _dropWidget->setShowDropIndicator(false);
        computeMetaData();
    }
}

const std::vector<unsigned int>& DEMultiSpeciesPlugin::getSpeciesIDs(const size_t species)
{
    _mappedSpeciesIDs.clear();

    if (_useSelMapForDE == 2) { // Map from species to points
        const auto otherData = _additionalSettingsDialog.getSelectionMappingSourcePicker().getCurrentDataset<Points>();
        const auto [selectionMapping, numPointsTarget] = getSelectionMappingOtherToCurrent(otherData, _points);
        const bool validSelectionMap = isSurjectiveMappingValid(selectionMapping, numPointsTarget, _points);

        if (!validSelectionMap)
            qDebug() << "ClusterDEMultiSpeciesPlugin: Invalid selectin map - things are about to break";

        const std::map<std::uint32_t, std::vector<std::uint32_t>>& mapSpeciesToPoints = selectionMapping->getMapping().getMap();

        const auto& speciesIndices = _clusters->getClusters()[species].getIndices();
        for (const std::uint32_t speciesID : speciesIndices) {

            if (!mapSpeciesToPoints.contains(speciesID))
                continue;

            const auto& mappedIDs = mapSpeciesToPoints.at(speciesID);
            _mappedSpeciesIDs.insert(_mappedSpeciesIDs.end(), mappedIDs.begin(), mappedIDs.end());
        }

        local::sortAndUnique(_mappedSpeciesIDs);

        return _mappedSpeciesIDs;
    }

    auto& clusterIDs = _clusters->getClusters()[species].getIndices();
    local::sortAndUnique(clusterIDs);

    return clusterIDs;
}

void DEMultiSpeciesPlugin::computeMetaData()
{
    if (!(_points.isValid() && _clusters.isValid()))
        return;

    // Compute normalization
    const auto& speciesClusters = _clusters->getClusters();
    const auto numSpecies       = speciesClusters.size();
    const auto numDimensions    = _points->getNumDimensions();
    const auto numPoints        = _points->getNumPoints();

    qDebug() << "ClusterDEMultiSpeciesPlugin: Computing dimension ranges";
    local::resizeNestedVec(_minValues, numSpecies, numDimensions, std::numeric_limits<float>::max());
    local::resizeNestedVec(_rescaleValues, numSpecies, numDimensions, std::numeric_limits<float>::max());
    
    int species = 0;

    auto computeMinAndRescale = [this, &species](auto globalRowID, auto localRowID, auto column, auto value) -> void 
        {
            if (value > _rescaleValues[species][column])
                _rescaleValues[species][column] = value;

            if (value < _minValues[species][column])
                _minValues[species][column] = value;
        };

    for (; species < numSpecies; species++) {
        const std::vector<unsigned int>& speciesIDs = getSpeciesIDs(species);

        local::visitElements(_points, speciesIDs, computeMinAndRescale);

        // Compute rescale values
#pragma omp parallel for schedule(dynamic,1)
        for (std::ptrdiff_t dim = 0; dim < numDimensions; dim++)
        {
            const float diff = (_rescaleValues[species][dim] - _minValues[species][dim]);
            if (std::fabs(diff) > 1e-6f)
                _rescaleValues[species][dim] = 1.0f / diff;
            else
                _rescaleValues[species][dim] = 1.0f;
        }
    }

    qDebug() << "DEMultiSpeciesPlugin: Loaded " << numDimensions << " dimensions for " << numPoints << " points";
}

void DEMultiSpeciesPlugin::writeToCSV() const
{
    if (_tableItemModel.isNull())
        return;

    // Let the user chose the save path
    QSettings settings(QLatin1String{ "ManiVault" }, QLatin1String{ "Plugins/" } + getKind());
    const QLatin1String directoryPathKey("directoryPath");
    const auto directoryPath = settings.value(directoryPathKey).toString() + "/";

    QString fileName = QFileDialog::getSaveFileName(
        nullptr, tr("Save data set"), directoryPath + "DifferentialExpression.csv", tr("CSV file (*.csv);;All Files (*)"));

    // Only continue when the dialog has not been not canceled and the file name is non-empty.
    if (fileName.isNull() || fileName.isEmpty())
    {
        qDebug() << "ClusterDEMultiSpeciesPlugin: No data written to disk - File name empty";
        return;
    }
    else
    {
        // store the directory name
        settings.setValue(directoryPathKey, QFileInfo(fileName).absolutePath());
    }

    QString csvString = _tableItemModel->createCSVString(',');
    if (csvString.isEmpty())
        return;
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Truncate))
        return;
    QTextStream output(&file);
    output << csvString;
    file.close();
}

void DEMultiSpeciesPlugin::computeDE()
{
    if (!(_points.isValid() && _clusters.isValid()))
        return;

    _tableItemModel->invalidate();

    // Compute differential expr
    qDebug() << "ClusterDEMultiSpeciesPlugin: Computing differential expression.";

    auto& speciesClusters               = _clusters->getClusters();
    const size_t numSpecies             = speciesClusters.size();
    const std::ptrdiff_t numDimensions = _points->getNumDimensions();
    const size_t selectionSizeA        = _selectionA.size();
    const size_t selectionSizeB        = _selectionB.size();

    // for mean, sum all values and divide by size later
    std::vector<std::vector<float>> meansA;
    std::vector<std::vector<float>> meansB;

    local::resizeNestedVec(meansA, numSpecies, numDimensions, 0.f);
    local::resizeNestedVec(meansB, numSpecies, numDimensions, 0.f);

    auto computeAvg= [this](const std::vector<uint32_t>& selectionIDs, std::vector<float>& means) -> void {
        local::visitElements(_points, selectionIDs, [&means](auto globalRowID, auto localRowID, auto column, auto value) -> void {
                means[column] += value;
            });
        };

    auto normAvg = [&](const std::vector<float>& avgs, const std::vector<float>& mins, const std::vector<float>& norm, const std::ptrdiff_t dim) -> float {
        return (avgs[dim] - mins[dim]) * norm[dim];
        };

    auto intersection = [](const std::vector<unsigned int>& a, const std::vector<unsigned int>& b) -> std::vector<unsigned int> {
        std::vector<unsigned int> intersection;
        std::set_intersection(
            a.cbegin(), a.cend(),
            b.cbegin(), b.cend(),
            std::back_inserter(intersection)
        );
        return intersection;
        };

    for (size_t species = 0; species < numSpecies; species++) {
        const std::vector<unsigned int>& speciesID = getSpeciesIDs(species);

        auto& meansA_species = meansA[species];
        auto& meansB_species = meansB[species];

        // We know that _selectionA, _selectionB and speciesID are sorted, it's done above
        const std::vector<unsigned int> selA_species = intersection(_selectionA, speciesID);
        const std::vector<unsigned int> selB_species = intersection(_selectionB, speciesID);

        computeAvg(selA_species, meansA_species);
        computeAvg(selB_species, meansB_species);

        const auto& mins_species = _minValues[species];
        const auto& norms_species = _rescaleValues[species];

#pragma omp parallel for schedule(dynamic,1)
        for (std::ptrdiff_t d = 0; d < numDimensions; d++)
        {
            // first divide means by number of rows
            meansA_species[d] /= selectionSizeA;
            meansB_species[d] /= selectionSizeB;

            // then min max - optional by toggle action
            if (_norm) {
                meansA_species[d] = normAvg(meansA_species, mins_species, norms_species, d);
                meansB_species[d] = normAvg(meansB_species, mins_species, norms_species, d);
            }
        }

    }

    const auto& dimensionNames = _points->getDimensionNames();

    _tableItemModel->startModelBuilding(_totalTableColumns, numDimensions);
#pragma omp  parallel for schedule(dynamic,1)
    for (std::ptrdiff_t dimension = 0; dimension < numDimensions; ++dimension)
    {
        std::vector<QVariant> dataVector = { dimensionNames[dimension] };
        dataVector.reserve(_totalTableColumns);
            
        for (size_t species = 0; species < numSpecies; species++) {
            dataVector.push_back(local::fround(meansA[species][dimension] - meansB[species][dimension], 3));    // Differential expression
            dataVector.push_back(local::fround(meansA[species][dimension], 3));
            dataVector.push_back(local::fround(meansB[species][dimension], 3));
        }

        assert(dataVector.size() == _totalTableColumns);

        _tableItemModel->setRow(dimension, dataVector, Qt::Unchecked, true);
    }

    _tableItemModel->endModelBuilding();
}

void DEMultiSpeciesPlugin::tableView_clicked(const QModelIndex& index)
{
    if (_tableItemModel->status() != TableModel::Status::UpToDate)
        return;

    try
    {
        const QModelIndex firstColumn = index.sibling(index.row(), 0);
        const QString dimensionName = firstColumn.data().toString();

        _selectedIdAction.setString(dimensionName);
        _currentSelectedDimension.setCurrentDimensionName(dimensionName);
    }
    catch (...) // catch everything
    {
        qDebug() << "DEMultiSpeciesPlugin::tableView_clicked -> something went wrong :(";
    }
}

void DEMultiSpeciesPlugin::tableView_selectionChanged(const QItemSelection& selected, const QItemSelection& deselected)
{
    tableView_clicked(selected.indexes().first());
}

/******************************************************************************
 * Serialization
 ******************************************************************************/

void DEMultiSpeciesPlugin::fromVariantMap(const QVariantMap& variantMap)
{
    ViewPlugin::fromVariantMap(variantMap);

    for (auto action : _serializedActions)
    {
        if (variantMap.contains(action->getSerializationName()))
            action->fromParentVariantMap(variantMap);

    }

    _additionalSettingsDialog.fromParentVariantMap(variantMap);
    _clustersDataGUID.fromParentVariantMap(variantMap);

    QVariantMap propertiesMap = local::get_strict_value<QVariantMap>(variantMap.value("#Properties"));
    if (!propertiesMap.isEmpty())
    {
        {
            auto found = propertiesMap.constFind("TableViewHeaderState");
            if (found != propertiesMap.constEnd())
            {

                QVariant value = found.value();
                QString stateAsQString = local::get_strict_value<QString>(value);
                // When reading a QByteArray back from jsondocument it's a QString. to Convert it back to a QByteArray we need to use .toUtf8().
                QByteArray state = QByteArray::fromBase64(stateAsQString.toUtf8());
                assert(local::is_valid_QByteArray(state));
                _headerState = state;
            }
        }
    }

    setPositionDataset(_points, false);

    if (!_clustersDataGUID.getString().isEmpty()) {
        setClustersDataset(mv::data().getDataset(_clustersDataGUID.getString()));
    }

}

QVariantMap DEMultiSpeciesPlugin::toVariantMap() const
{
    QVariantMap variantMap = ViewPlugin::toVariantMap();

    for (auto action : _serializedActions)
    {
        assert(action->getSerializationName() != "#Properties");
        action->insertIntoVariantMap(variantMap);
    }

    _additionalSettingsDialog.insertIntoVariantMap(variantMap);
    _clustersDataGUID.insertIntoVariantMap(variantMap);

    // properties map
    QVariantMap propertiesMap;

    QByteArray headerState = _tableView->horizontalHeader()->saveState();
    propertiesMap["TableViewHeaderState"] = QString::fromUtf8(headerState.toBase64()); // encode the state with toBase64() and put it in a Utf8 QString since it will do that anyway. Best to be explicit in case it changes in the future
    variantMap["#Properties"] = propertiesMap;


    return variantMap;
}

// =============================================================================
// Factory
// =============================================================================

DEMultiSpeciesPluginFactory::DEMultiSpeciesPluginFactory()
{
    setIconByName("table");

    //for (const auto& tutorial_file : list_tutorial_files("tutorials/DEMultiSpecies")) {
    //    if (insert_md_into_json(tutorial_file)) {

    //        if (auto tutorial_json = readJSON(tutorial_file)) {
    //            mv::help().addTutorial(new LearningCenterTutorial(tutorial_json.value()["tutorials"].toArray().first().toObject().toVariantMap()));
    //        }
    //        
    //    }
    //}
}

ViewPlugin* DEMultiSpeciesPluginFactory::produce()
{
    return new DEMultiSpeciesPlugin(this);
}

mv::DataTypes DEMultiSpeciesPluginFactory::supportedDataTypes() const
{
    return { PointType } ;
}

mv::gui::PluginTriggerActions DEMultiSpeciesPluginFactory::getPluginTriggerActions(const mv::Datasets& datasets) const
{
    PluginTriggerActions pluginTriggerActions;

    const auto getPluginInstance = [this]() -> DEMultiSpeciesPlugin* {
        return dynamic_cast<DEMultiSpeciesPlugin*>(plugins().requestViewPlugin(getKind()));
        };

    const auto numberOfDatasets = datasets.count();

    if (numberOfDatasets >= 1 && PluginFactory::areAllDatasetsOfTheSameType(datasets, PointType)) {
        auto pluginTriggerAction = new PluginTriggerAction(const_cast<DEMultiSpeciesPluginFactory*>(this), this, "DE (Multiple species)", "Compute differential expressions between two selections for multiple species", StyledIcon(), [this, getPluginInstance, datasets](PluginTriggerAction& pluginTriggerAction) -> void {
            for (const auto& dataset : datasets)
                getPluginInstance()->setPositionDataset( dataset );
            });

        pluginTriggerActions << pluginTriggerAction;
    }

    return pluginTriggerActions;
}
