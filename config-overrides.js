const webpack = require('webpack');

module.exports = function override(config) {
  const fallback = config.resolve.fallback || {};
  config.module.rules.unshift({
    test: /\.m?js$/,
    resolve: {
      fullySpecified: false, // disable the behavior
    },
  });
  Object.assign(fallback, {
    "buffer": require.resolve("buffer"),
    "stream": require.resolve("stream-browserify"),
    "process": require.resolve("process/browser"),
    "assert": require.resolve("assert/")
  })
  config.resolve.fallback = fallback;
  config.plugins = (config.plugins || []).concat([
    new webpack.ProvidePlugin({
      Buffer: ['buffer', 'Buffer'],
      process: 'process/browser',
    })
  ])
  return config;
}
