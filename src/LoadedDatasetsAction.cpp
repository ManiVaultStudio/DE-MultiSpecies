#include "LoadedDatasetsAction.h"
#include "DE-MultiSpeciesPlugin.h"

#include "PointData/PointData.h"

#include <QMenu>

using namespace mv;
using namespace mv::gui;

LoadedDatasetsAction::LoadedDatasetsAction(QObject* parent, const QString& title) :
    GroupAction(parent, "Loaded datasets"),
    _positionDatasetPickerAction(this, "Position")
{
    setIcon(mv::util::StyledIcon("database"));
    setToolTip("Manage loaded datasets for position");

    _positionDatasetPickerAction.setFilterFunction([](const mv::Dataset<DatasetImpl>& dataset) -> bool {
        if(dataset->getDataType() == PointType)
        {
            return true;
        }
        return false;
    });
}

void LoadedDatasetsAction::initialize(DEMultiSpeciesPlugin* plugin)
{
    Q_ASSERT(plugin != nullptr);

    if (plugin == nullptr)
        return;

    _plugin = plugin;

    connect(&_positionDatasetPickerAction, &DatasetPickerAction::datasetPicked, [this](Dataset<DatasetImpl> pickedDataset) -> void {
        _plugin->getPositionDataset() = pickedDataset;
    });

    connect(&_plugin->getPositionDataset(), &Dataset<Points>::changed, this, [this](DatasetImpl* dataset) -> void {
        _positionDatasetPickerAction.setCurrentDataset(dataset);
    });
}

void LoadedDatasetsAction::fromVariantMap(const QVariantMap& variantMap)
{
    WidgetAction::fromVariantMap(variantMap);

    _positionDatasetPickerAction.fromParentVariantMap(variantMap);

    // Load position dataset
    auto positionDataset = _positionDatasetPickerAction.getCurrentDataset();
    if (positionDataset.isValid())
    {
        Dataset pickedDataset = mv::data().getDataset(positionDataset.getDatasetId());
        _plugin->getPositionDataset() = pickedDataset;
    }
}

QVariantMap LoadedDatasetsAction::toVariantMap() const
{
    QVariantMap variantMap = GroupAction::toVariantMap();

    _positionDatasetPickerAction.insertIntoVariantMap(variantMap);

    return variantMap;
}
