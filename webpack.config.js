const path = require('path');
const webpack = require('webpack');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const { PyodidePlugin } = require("@pyodide/webpack-plugin");
const CopyWebpackPlugin = require('copy-webpack-plugin');
const Dotenv = require('dotenv-webpack');

module.exports = {
  entry: './src/index.js',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'bundle.js',
    publicPath: '/'
  },
  mode: 'development',
  target: 'electron-renderer',
  module: {
    rules: [
      {
        test: /\.(js|jsx)$/,
        exclude: /node_modules/,
        use: 'babel-loader'
      },
      {
        test: /\.css$/,
        use: [
          'style-loader',
          'css-loader',
          'postcss-loader',
        ],
      },
      {
        test: /\.(woff|woff2|eot|ttf|otf)$/i,
        type: 'asset/resource',
        generator: {
          filename: 'assets/fonts/[name][ext]'
        }
      },
      {
        test: /\.(png|svg|jpg|jpeg|gif)$/i, // Rule to handle image files
        use: [
          {
            loader: 'file-loader',
            options: {
              name: '[name].[ext]',
              outputPath: 'assets/', // Unified path for all assets
              publicPath: 'assets/' // Unified public path
            }
          }
        ]
      },
      {
        test: /\.worker\.js$/,
        use: {
          loader: 'worker-loader',
          options: {
            filename: '[name].[contenthash].worker.js'
          }
        }
      }
    ]
  },
  resolve: {
    extensions: [".js", ".jsx"],
    alias: {
      global: path.resolve(__dirname, "./global-shim.js"),
    },
    fallback: {
      "buffer": require.resolve("buffer/"),
      "stream": require.resolve("stream-browserify"),
      "process/browser": require.resolve("process/browser"),
      "path": require.resolve("path-browserify"),
      "crypto": require.resolve("crypto-browserify"),
    }
  },
  plugins: [
    new webpack.ProvidePlugin({
      Buffer: ['buffer', 'Buffer'],
      process: 'process/browser'
    }),
    new HtmlWebpackPlugin({
      template: './public/index.html',
      filename: 'index.html'
    }),
    new webpack.HotModuleReplacementPlugin(),
    new PyodidePlugin(),
    new CopyWebpackPlugin({
      patterns: [
        { from: 'public/py-wheels', to: 'py-wheels' }
        //{ from: 'public/wasm', to: 'wasm' }
      ],
    }),
    new Dotenv(),
  ],
  devServer: {
    static: {
      directory: path.join(__dirname, 'public')
    },
    compress: true,
    port: 9328,
    historyApiFallback: true,
    proxy: {
      '/datasets': {
        target: 'https://docgl1or94tw4.cloudfront.net',
        changeOrigin: true,
        pathRewrite: { '^/datasets': '' },
      },
      '/broadinstitute': {
        target: 'https://data.broadinstitute.org',
        changeOrigin: true,
        pathRewrite: { '^/broadinstitute': '' },
      },
      '/msigdb': {
        target: 'https://www.gsea-msigdb.org',
        changeOrigin: true,
        pathRewrite: { '^/msigdb': '' },
      },
    },
  }
};
