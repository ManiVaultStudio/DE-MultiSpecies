#pragma once

#include <actions/GroupAction.h>

#include "actions/DatasetPickerAction.h"

using namespace mv::gui;

class DEMultiSpeciesPlugin;

class LoadedDatasetsAction : public GroupAction
{
    Q_OBJECT

public:
    Q_INVOKABLE LoadedDatasetsAction(QObject* parent, const QString& title);

    void initialize(DEMultiSpeciesPlugin* plugin);

public: // Serialization

    /**
     * Load widget action from variant map
     * @param Variant map representation of the widget action
     */
    void fromVariantMap(const QVariantMap& variantMap) override;

    /**
     * Save widget action to variant map
     * @return Variant map representation of the widget action
     */
    QVariantMap toVariantMap() const override;

private:
    DEMultiSpeciesPlugin*   _plugin = nullptr;
    DatasetPickerAction     _positionDatasetPickerAction;

    friend class mv::AbstractActionsManager;
};

Q_DECLARE_METATYPE(LoadedDatasetsAction)

inline const auto loadedDatasetsActionMetaTypeId = qRegisterMetaType<LoadedDatasetsAction*>("LoadedDatasetsAction");
