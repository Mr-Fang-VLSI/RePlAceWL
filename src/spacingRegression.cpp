#include <shark/Algorithms/DirectSearch/CMA.h>
#include <shark/Core/Types.h>
#include <shark/Models/Trees/GeneralizedBoostedTrees.h>
#include <shark/ObjectiveFunctions/Loss/SquaredLoss.h>
#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/Algorithms/Trainers/GBTTrainer.h>
#include <shark/Data/Dataset.h>
#include <iostream>
#include <vector>

using namespace shark;

// 自定义优化目标函数，计算在给定spacing值时的预测wirelength
class WirelengthObjectiveFunction : public SingleObjectiveFunction {
public:
    WirelengthObjectiveFunction(const GBT<RealVector> &model, double sizeOfNetlist)
        : m_model(model), m_sizeOfNetlist(sizeOfNetlist) {
        m_features.resize(2);
        m_features[1] = sizeOfNetlist;
    }

    std::size_t numberOfVariables() const {
        return 1;  // 只优化一个变量，即spacing
    }

    double eval(const RealVector &spacing) const {
        m_features[0] = spacing[0];
        RealVector predictedWirelength = m_model(m_features);
        return predictedWirelength[0];
    }

private:
    const GBT<RealVector> &m_model;
    RealVector m_features;
    double m_sizeOfNetlist;
};

int main() {
    //... (之前的代码)
    std::vector<RealVector> features = {
        {3, 10},
        {2, 15},
        {4, 20},
        {1, 25},
        {5, 30},
    };
    std::vector<RealVector> wirelength_labels = {
        {100},
        {200},
        {300},
        {400},
        {500},
    };

    // 创建数据集
    Data<RealVector> inputs = createDataFromRange(features);
    Data<RealVector> labels = createDataFromRange(wirelength_labels);
    ClassificationDataset dataset = createLabeledDataFromRange(inputs, labels);

    // 划分训练集和测试集
    ClassificationDataset trainSet, testSet;
    tie(trainSet, testSet) = splitAtElement(dataset, 4);

    // 创建梯度提升树模型
    GBTTrainer<RealVector> trainer;
    trainer.setNumberOfBoostingIterations(100);
    trainer.setBaseModelCount(100);
    trainer.setRegularization(0.1);
    SquaredLoss<RealVector> loss;

    // 训练模型
    GBT<RealVector> model;
    trainer.train(model, trainSet, loss);

    // 使用模型预测测试数据
    std::size_t correctPredictions = 0;
    for (std::size_t i = 0; i < testSet.numberOfElements(); ++i) {
        RealVector prediction = model(testSet.element(i).input);
        std::cout << "Prediction: " << prediction << " - Actual: " << testSet.element(i).label << std::endl;
    }

    // 预测新数据
    std::vector<RealVector> new_features = {
        {2, 18},  // spacing = 2, size_of_netlist = 18
        {4, 12},  // spacing = 4, size_of_netlist = 12
    };

    for (const auto &feature : new_features) {
        RealVector predicted_wirelength = model(feature);
        std::cout << "Predicted wirelength: " << predicted_wirelength << std::endl;
    }

    // 使用 CMA-ES 优化器找到最佳的 spacing 值
    WirelengthObjectiveFunction objectiveFunction(model, 16);
    CMA optimizer;

    // 初始化优化器
    RealVector initialSpacing(1, 3.0); // 可以设置一个合理的初始值
    optimizer.init(objectiveFunction, initialSpacing);

    // 运行优化
    for (std::size_t step = 0; step < 100; ++step) {
        optimizer.step(objectiveFunction);
    }

    // 输出最优解
    double bestSpacing = optimizer.solution().point[0];
    std::cout << "Best spacing for size_of_netlist = 16: " << bestSpacing << std::endl;
    std::cout << "Minimum wirelength: " << optimizer.solution().value << std::endl;

    return 0;
}
