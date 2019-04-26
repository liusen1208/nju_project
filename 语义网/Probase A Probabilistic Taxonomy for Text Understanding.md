<p align="center">
    <h1 align="center">Probase: A Probabilistic Taxonomy for Text Understanding
</h1>
    <p align="center">A universal, probabilistic taxonomy which uses probabilities to model inconsistent, ambiguous and uncertain information it contains.<br><br>一种使用概率来模拟其包含的不一致,模糊和不确定信息的概念分类法</p>
<br>
</p>

## Probase的介绍

1. Probase 三个独特之处：
   * 采用了新框架，包含了一个迭代学习算法和一个分类构造算法组成。
   * 使用一种概率的方法来模拟它拥有的知识，每个事实或关系都与一些概率相关联，以衡量其合理性和典型性。这种概率处理使Probase能更好的捕获人类语言的语意。
   * 它完全由Web上的HTML文本自动构建，是最大的通用分类法。其概念空间巨大，比YAGO大8倍。

## 迭代提取

知识获取包含信息提取，数据处理和整合。 在信息提取过程中，Probase并没有选择常用的句法迭代(Syntactic Iteration)， 而是结合了句法模式和语义模式创造出的一种新的方式来提取信息。

### 新的框架

为了针对单句法模式(Syntactic Patterns)的不能准确提取名词短语的不足，文章提出一种新的框架

1. 句法提取：
   * 句法提取：

